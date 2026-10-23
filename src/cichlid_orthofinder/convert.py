"""OrthoFinder's orthogroups and gene trees -> the project's assign/ and families/ (sheet sections 5 and 6).

    cichlid-orthofinder convert --results <OrthoFinder Results_<date>/> --run-settings <run_settings.json>
        [--tree-set raw_pruned|resolved] [--alignment orthogroup|raw_cluster] [--min-size 4]
        [--cache <tree_cache>] [--out <Results/v1/OrthoFinder/>]

The converter reads explicit input and output paths, so it can be tested on a scratch tree; by default the
command line takes the output folder from the coordinator (`get_layout().results(release, "OrthoFinder")`).
Every format is written with the coordinator's writers (`schemas.ASSIGN.write`, `schemas.write_run`) and
checked with `cichlid-data validate`.

Choices that are the owner's, so they are parameters here (defaults as recommended, not decided):

- tree set: `raw_pruned`, the raw IQ-TREE tree of the preliminary cluster (`WorkingDirectory/Trees_ids/`)
  pruned to the orthogroup's genes, unrooted (the coordinator's recommendation); or `resolved`, OrthoFinder's
  rooted `Resolved_Gene_Trees/Resolved_Gene_Trees.txt`;
- alignment source: `orthogroup`, `MultipleSequenceAlignments/<OG>.fa` (the sheet's wording); or
  `raw_cluster`, the preliminary cluster's alignment (`Alignments_ids/`, the one the raw tree was built on),
  restricted to the orthogroup's rows with the columns that are then all gaps removed;
- minimum family size (genes): 4.

Determinism: the same OrthoFinder output gives the same files, byte for byte. Rows and families are
sorted; branch lengths keep IQ-TREE's strings where a length is unchanged and sums are written with
"%.10g"; lines end in "\\n"; family files hold no timestamps.
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

TREE_SETS = ("raw_pruned", "resolved")
ALIGNMENT_SOURCES = ("orthogroup", "raw_cluster")
DEFAULT_MIN_SIZE = 4
RELEASE, PART, METHOD = "v1", "OrthoFinder", "OrthoFinder"
INTERNAL_RE = re.compile(r"^\d+_\d+$")
ROOTING = "orthofinder species tree"


class ConversionError(ValueError):
    """The OrthoFinder output cannot be converted as it is: the message says why."""


# ------------------------------------------------------------------------------------------------- Newick

class Node:
    """A tree node. `length` is the branch length as written (a string), or None."""
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name="", length=None):
        self.name, self.length, self.children, self.parent = name, length, [], None

    def add(self, child):
        child.parent = self
        self.children.append(child)
        return child

    def iter(self):
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def leaves(self):
        return [n for n in self.iter() if not n.children]


def parse_newick(text):
    """Parse a Newick tree, keeping each branch length's text. Iterative (deep trees)."""
    s = text.strip()
    if not s.endswith(";"):
        raise ConversionError("a Newick tree must end with ';'")
    s = s[:-1]
    root = cur = Node()
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            cur = cur.add(Node())
            i += 1
        elif c == ",":
            if cur.parent is None:
                raise ConversionError("unbalanced parentheses in a tree")
            cur = cur.parent.add(Node())
            i += 1
        elif c == ")":
            if cur.parent is None:
                raise ConversionError("unbalanced parentheses in a tree")
            cur = cur.parent
            i += 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in "(),:;":
                j += 1
            cur.length = s[i + 1:j].strip() or None
            i = j
        elif c.isspace():
            i += 1
        else:
            j = i
            while j < n and s[j] not in "(),:;":
                j += 1
            cur.name = s[i:j].strip()
            i = j
    if cur is not root:
        raise ConversionError("unbalanced parentheses in a tree")
    return root


def write_newick(root):
    """Leaves with their names, internal nodes without names (OrthoFinder's n0, n1, ... are dropped); the
    root has no branch length. Iterative."""
    out, stack = [], [("open", root)]
    while stack:
        kind, node = stack.pop()
        if kind == "text":
            out.append(node)
        elif node.children:
            out.append("(")
            stack.append(("text", ")" + (f":{node.length}" if node.length is not None and node is not root else "")))
            for k, child in enumerate(reversed(node.children)):
                stack.append(("open", child))
                if k < len(node.children) - 1:
                    stack.append(("text", ","))
        else:
            out.append(node.name + (f":{node.length}" if node.length is not None and node is not root else ""))
    return "".join(out) + ";"


def _add_lengths(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return "%.10g" % (float(a) + float(b))


def _detach(node):
    node.parent.children.remove(node)
    node.parent = None


def prune(root, keep):
    """Keep only the leaves named in `keep`; remove the nodes left with one child, adding the branch lengths
    of the two branches they joined. Returns the (possibly new) root."""
    for leaf in root.leaves():
        if leaf.name not in keep and leaf is not root:
            _detach(leaf)
    changed = True
    while changed:                                     # internal nodes left without children
        changed = False
        for node in list(root.iter()):
            if node is not root and not node.children and node.name not in keep:
                _detach(node)
                changed = True
    for node in list(root.iter()):                     # degree-2 nodes
        if node is not root and len(node.children) == 1:
            child, parent = node.children[0], node.parent
            child.length = _add_lengths(child.length, node.length)
            parent.children[parent.children.index(node)] = child
            child.parent = parent
    while len(root.children) == 1:
        root = root.children[0]
        root.parent, root.length = None, None
    return root


def unroot(root):
    """A root with two children becomes a root with three (the unrooted form): the two root branches are
    merged into one, their lengths added. A two-leaf tree is left as it is."""
    if len(root.children) != 2:
        return root
    a, b = root.children
    inner, other = (a, b) if a.children else (b, a)
    if not inner.children:
        return root
    other.length = _add_lengths(other.length, inner.length)
    inner.parent, inner.length = None, None
    root.children = []
    inner.add(other)
    return inner


def path_lengths(root):
    """-> {(leaf a, leaf b): summed branch length between them} for every pair (for tests and checks)."""
    depth, up = {}, {}
    for node in root.iter():
        depth[node] = 0.0 if node.parent is None else depth[node.parent] + float(node.length or 0)
    leaves = root.leaves()
    ancestors = {}
    for leaf in leaves:
        chain, n = [], leaf
        while n is not None:
            chain.append(n)
            n = n.parent
        ancestors[leaf] = chain
    out = {}
    for i, a in enumerate(leaves):
        sa = set(ancestors[a])
        for b in leaves[i + 1:]:
            lca = next(n for n in ancestors[b] if n in sa)
            out[tuple(sorted((a.name, b.name)))] = depth[a] + depth[b] - 2 * depth[lca]
    return out


# ------------------------------------------------------------------------------------------------- inputs

def read_fasta(path):
    """-> [(name, sequence)] in file order (the name is the header's first word)."""
    rows, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, "".join(seq)))
                name, seq = (line[1:].split() or [""])[0], []
            elif name is not None:
                seq.append(line.strip())
    if name is not None:
        rows.append((name, "".join(seq)))
    return rows


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_species_ids(working):
    """SpeciesIDs.txt -> {species number: prefix}. OrthoFinder numbers species alphabetically: read, never assume."""
    out = {}
    for line in Path(working, "SpeciesIDs.txt").read_text().splitlines():
        if line.strip():
            number, name = line.split(":", 1)
            out[int(number)] = name.strip().removesuffix(".faa")
    return out


def read_sequence_ids(working):
    """SequenceIDs.txt -> {internal name (0_393): our gene ID}."""
    out = {}
    for line in Path(working, "SequenceIDs.txt").read_text().splitlines():
        if line.strip():
            internal, name = line.split(":", 1)
            out[internal.strip()] = name.strip().split()[0]
    return out


def clusters_file(working):
    """The preliminary MCL clusters with OrthoFinder's internal names (clusters_OrthoFinder_I*.txt_id_pairs.txt)."""
    found = sorted(Path(working).glob("clusters_OrthoFinder_I*.txt_id_pairs.txt"))
    if len(found) != 1:
        raise ConversionError(f"expected one clusters_OrthoFinder_I*.txt_id_pairs.txt in {working}, found {len(found)}")
    return found[0]


def read_clusters(working):
    """-> {cluster number: [internal names]}. Cluster k's tree is Trees_ids/OG<k, 7 digits>.txt and its
    alignment Alignments_ids/OG<k>.fa (checked on pilot 2: 512 of 512)."""
    text = clusters_file(working).read_text()
    if "begin" not in text:
        raise ConversionError("the clusters file has no MCL matrix")
    body = text.split("begin", 1)[1].rsplit(")", 1)[0]
    out = {}
    for chunk in body.split("$"):
        tokens = chunk.split()
        if tokens:
            out[int(tokens[0])] = tokens[1:]
    return out


def read_orthogroups(results):
    """Orthogroups/Orthogroups.txt -> {orthogroup: [gene IDs]}, and the single-gene orthogroups of
    Orthogroups_UnassignedGenes.tsv as a set of genes (a gene there is in no orthogroup)."""
    groups = {}
    for line in Path(results, "Orthogroups", "Orthogroups.txt").read_text().splitlines():
        if line.strip():
            og, genes = line.split(":", 1)
            groups[og.strip()] = sorted(genes.split())
    unassigned = set()
    with open(Path(results, "Orthogroups", "Orthogroups_UnassignedGenes.tsv")) as fh:
        next(fh, None)
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            unassigned.update(c for c in cells[1:] if c)
    return {og: g for og, g in groups.items() if not (len(g) == 1 and g[0] in unassigned)}, unassigned


def strip_species(label, prefixes):
    """Acal_Acal_113012997 -> Acal_113012997: exactly one repeated <prefix>_ is removed."""
    for p in prefixes:
        if label.startswith(f"{p}_{p}_"):
            return label[len(p) + 1:]
    raise ConversionError(f"label {label!r} does not start with a doubled species prefix")


# ------------------------------------------------------------------------------------------------- outputs

class OutputPaths:
    """Where the conversion writes. For the real results, the coordinator's ResultsLayout
    (`get_layout().results(release, "OrthoFinder")`): every path comes from its methods. For tests and samples,
    the same structure under a scratch folder."""

    def __init__(self, results_layout=None, root=None):
        if (results_layout is None) == (root is None):
            raise ValueError("give a ResultsLayout or a scratch folder, not both")
        self.layout, self.root = results_layout, Path(root) if root is not None else None

    @classmethod
    def of(cls, out):
        """A ResultsLayout, an OutputPaths, or a scratch folder -> OutputPaths."""
        if isinstance(out, cls):
            return out
        if hasattr(out, "family") and hasattr(out, "run_record"):
            return cls(results_layout=out)
        return cls(root=out)

    def assign(self, prefix):
        return Path(self.layout.assign(prefix).local) if self.layout else self.root / "assign" / f"{prefix}.assign.tsv"

    def families(self):
        return Path(self.layout.section("families").local) if self.layout else self.root / "families"

    def family(self, og):
        return Path(self.layout.family(og).local) if self.layout else self.root / "families" / og

    def native(self):
        return Path(self.layout.native().local) if self.layout else self.root / "native"

    def run_record(self):
        return Path(self.layout.run_record().local) if self.layout else self.root / "run.json"


# ------------------------------------------------------------------------------------------------- one family

def raw_tree(working, cluster, genes, internal_to_gene):
    """The cluster's raw IQ-TREE tree, leaves renamed to our IDs, pruned to `genes`, unrooted."""
    f = Path(working, "Trees_ids", f"OG{cluster:07d}.txt")
    if not f.is_file():
        return None
    root = parse_newick(f.read_text())
    for leaf in root.leaves():
        if leaf.name not in internal_to_gene:
            raise ConversionError(f"{f.name}: leaf {leaf.name!r} is not in SequenceIDs.txt")
        leaf.name = internal_to_gene[leaf.name]
    root = unroot(prune(root, set(genes)))
    if {l.name for l in root.leaves()} != set(genes):
        raise ConversionError(f"{f.name}: the tree does not hold every gene of the orthogroup")
    return root


def resolved_trees(results):
    """Resolved_Gene_Trees.txt -> {orthogroup: tree text}."""
    out = {}
    f = Path(results, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    if f.is_file():
        for line in f.read_text().splitlines():
            if line.strip():
                og, tree = line.split(":", 1)
                out[og.strip()] = tree.strip()
    return out


def resolved_tree(text, genes, prefixes):
    root = parse_newick(text)
    for leaf in root.leaves():
        leaf.name = strip_species(leaf.name, prefixes)
    if sorted(l.name for l in root.leaves()) != sorted(genes):
        raise ConversionError("the resolved tree's leaves are not the orthogroup's genes")
    return root


def check_alignment(rows, leaves, where):
    """Refuse an alignment with a missing row for a leaf, a duplicate row or unequal row lengths."""
    names = [n for n, _ in rows]
    dup = sorted(n for n, c in Counter(names).items() if c > 1)
    if dup:
        raise ConversionError(f"{where}: duplicate rows {', '.join(dup[:5])}")
    missing = sorted(set(leaves) - set(names))
    if missing:
        raise ConversionError(f"{where}: no row for leaves {', '.join(missing[:5])}"
                              + (f" (rows: {', '.join(names[:3])})" if len(names) <= 3 else ""))
    lengths = {len(s) for n, s in rows if n in leaves}
    if len(lengths) != 1:
        raise ConversionError(f"{where}: rows of different lengths ({min(lengths)} to {max(lengths)})")


def alignment(results, source, og, cluster, genes, prefixes, gene_to_internal):
    """-> [(gene ID, aligned sequence)], sorted by gene ID, exactly one row per gene."""
    if source == "orthogroup":
        f = Path(results, "MultipleSequenceAlignments", f"{og}.fa")
        if not f.is_file():
            raise ConversionError(f"{og}: no MultipleSequenceAlignments/{og}.fa")
        rows = []
        for name, seq in read_fasta(f):
            try:
                rows.append((strip_species(name, prefixes), seq))
            except ConversionError:
                raise ConversionError(f"MultipleSequenceAlignments/{og}.fa: row {name!r} is not OrthoFinder's "
                                      f"<prefix>_<prefix>_<id> (rows: {len(genes)} genes expected)") from None
        check_alignment(rows, genes, f"MultipleSequenceAlignments/{og}.fa")
    else:
        f = Path(results, "WorkingDirectory", "Alignments_ids", f"OG{cluster:07d}.fa")
        if not f.is_file():
            raise ConversionError(f"{og}: no Alignments_ids/OG{cluster:07d}.fa")
        internal = {gene_to_internal[g]: g for g in genes}
        raw = read_fasta(f)
        check_alignment([(internal.get(n, n), s) for n, s in raw], genes, f"Alignments_ids/{f.name}")
        rows = [(internal[n], s) for n, s in raw if n in internal]
        if not rows:
            raise ConversionError(f"Alignments_ids/{f.name}: no rows for {og}")
        keep = [i for i in range(len(rows[0][1])) if any(s[i] != "-" for _, s in rows)]
        rows = [(n, "".join(s[i] for i in keep)) for n, s in rows]
    by_gene = dict(rows)
    return [(g, by_gene[g]) for g in sorted(genes)]


def inputs_hash(tree_text, aln_bytes, of_version, mapping_md5):
    """SHA-256 over the tree text, the alignment bytes, the OrthoFinder version and the mapping files' md5."""
    h = hashlib.sha256()
    for part in (tree_text.encode(), aln_bytes, of_version.encode(), mapping_md5.encode()):
        h.update(part)
        h.update(b"\n")
    return h.hexdigest()


def write_family(final, tree_root, rows, release, of_version, mapping_md5):
    """Write the family folder `final`: tree.nwk, combined.aln, then done.json last. -> "built", "reused" or
    "rebuilt". A finished family (done.json present) whose done.json matches, inputs hash included, and whose
    tree.nwk and combined.aln are byte for byte what would be written, is reused as it is; a folder without done.json is unfinished and is built again; a family whose inputs
    changed is rebuilt. Each is built in a hidden folder (which validation skips) and renamed into place, so a
    half-written family is never visible."""
    final = Path(final)
    tree_text = write_newick(tree_root) + "\n"
    aln = "".join(f">{g}\n{s}\n" for g, s in rows).encode()
    rooted = len(tree_root.children) == 2
    done = {"release": release, "rooted": rooted, "tree_mode": "denovo",
            "inputs_hash": inputs_hash(tree_text, aln, of_version, mapping_md5)}
    if rooted:
        done["rooting"] = ROOTING
    how = "built"
    if final.exists():
        try:
            old = json.loads((final / "done.json").read_text())
        except (OSError, ValueError):
            old = None                                   # unfinished: no done.json, or an unreadable one
        same_files = ((final / "tree.nwk").is_file() and (final / "tree.nwk").read_text() == tree_text
                      and (final / "combined.aln").is_file() and (final / "combined.aln").read_bytes() == aln)
        if old == done and same_files:
            return "reused"
        how = "rebuilt"
    tmp, aside = final.with_name(f".{final.name}.part"), final.with_name(f".{final.name}.old")
    for leftover in (tmp, aside):
        if leftover.exists():
            shutil.rmtree(leftover)
    tmp.mkdir(parents=True)
    (tmp / "tree.nwk").write_text(tree_text)
    (tmp / "combined.aln").write_bytes(aln)
    (tmp / "done.json").write_text(json.dumps(done, sort_keys=True, indent=1) + "\n")
    if final.exists():
        final.rename(aside)
    tmp.rename(final)
    if aside.exists():
        shutil.rmtree(aside)
    return how


# ------------------------------------------------------------------------------------------------- run

def models_by_cluster(cache):
    """-> {cluster ID (OG0000123): (sequences, model)} from the tree wrapper's events (the last 'done')."""
    out = {}
    events = Path(cache, "events.jsonl") if cache else None
    if events and events.is_file():
        for line in events.read_text().splitlines():
            try:
                e = json.loads(line)
            except ValueError:
                continue
            if e.get("event") == "done":
                out[e["og"]] = (e.get("sequences"), e.get("model", "MFP"))
    return out


def convert(results, out, species, proteomes, of_version, tree_set="raw_pruned", alignment_source="orthogroup",
            min_size=DEFAULT_MIN_SIZE, release=RELEASE, cache=None, say=print):
    """Convert one OrthoFinder Results_<date>/ folder into `out` (assign/, families/, native/summary/): the
    coordinator's ResultsLayout for the real results, or a scratch folder with the same structure. `species`
    is the catalog's ordered list of prefixes and `proteomes` {prefix: path}. A rerun into the same output
    reuses the finished families whose inputs are unchanged. -> a summary dict (counts and the split-cluster
    check) for run.json."""
    from cichlid_data import schemas
    if tree_set not in TREE_SETS:
        raise ConversionError(f"tree set must be one of {TREE_SETS}")
    if alignment_source not in ALIGNMENT_SOURCES:
        raise ConversionError(f"alignment source must be one of {ALIGNMENT_SOURCES}")
    results, paths = Path(results), OutputPaths.of(out)
    working = results / "WorkingDirectory"

    # 1. mapping files: every species and every sequence accounted for
    species_ids = read_species_ids(working)
    if sorted(species_ids.values()) != sorted(species):
        raise ConversionError(f"SpeciesIDs.txt has {sorted(species_ids.values())}, the catalog {sorted(species)}")
    internal_to_gene = read_sequence_ids(working)
    gene_to_internal = {g: i for i, g in internal_to_gene.items()}
    if len(gene_to_internal) != len(internal_to_gene):
        raise ConversionError("SequenceIDs.txt maps two sequences to the same gene ID")
    proteome_rows = {p: [(n, len(s)) for n, s in read_fasta(proteomes[p])] for p in species}
    n_proteins = sum(len(r) for r in proteome_rows.values())
    if n_proteins != len(internal_to_gene):
        raise ConversionError(f"SequenceIDs.txt has {len(internal_to_gene)} sequences, the proteomes {n_proteins}")
    for number, prefix in species_ids.items():
        theirs = {g for i, g in internal_to_gene.items() if i.split("_")[0] == str(number)}
        mine = {n for n, _ in proteome_rows[prefix]}
        if theirs != mine:
            raise ConversionError(f"{prefix}: SequenceIDs.txt and the proteome hold different genes "
                                  f"({len(theirs - mine)} only in OrthoFinder, {len(mine - theirs)} only in the proteome)")
    mapping_md5 = md5_file(working / "SpeciesIDs.txt") + md5_file(working / "SequenceIDs.txt")
    prefixes = sorted(species, key=len, reverse=True)       # longest first, for the doubled-prefix strip

    # 2. orthogroups, and the cluster that holds each
    groups, unassigned = read_orthogroups(results)
    clusters = read_clusters(working)
    cluster_of_internal = {i: k for k, members in clusters.items() for i in members}
    cluster_of = {}
    for og, genes in groups.items():
        found = {cluster_of_internal.get(gene_to_internal.get(g)) for g in genes}
        if len(found) != 1 or None in found:
            raise ConversionError(f"{og} does not lie in exactly one preliminary cluster ({sorted(map(str, found))})")
        cluster_of[og] = found.pop()
    in_group = {g: og for og, genes in groups.items() for g in genes}
    stray = sorted(set(gene_to_internal) - set(in_group) - unassigned)
    if stray:
        raise ConversionError(f"{len(stray)} genes are neither in an orthogroup nor unassigned, e.g. {stray[:3]}")

    # 3. assign/: every gene of each proteome, in the proteome's order
    for prefix in species:
        rows = [{"gene_id": g, "protein_length": length, "best_family": in_group.get(g, "")}
                for g, length in proteome_rows[prefix]]
        schemas.ASSIGN.write(paths.assign(prefix), rows, prefix)

    # 4. families/: which orthogroups get one; a family folder this conversion would not write is refused
    resolved = resolved_trees(results) if tree_set == "resolved" else {}
    small = sorted(og for og in groups if len(groups[og]) < min_size)
    has_tree = (lambda og: Path(working, "Trees_ids", f"OG{cluster_of[og]:07d}.txt").is_file()) \
        if tree_set == "raw_pruned" else (lambda og: og in resolved)
    no_tree = sorted(og for og in groups if og not in small and not has_tree(og))
    families = paths.families()
    expected = sorted(og for og in groups if og not in small and og not in no_tree)
    present = sorted(e.name for e in families.iterdir() if e.is_dir() and not e.name.startswith(".")) \
        if families.is_dir() else []
    stale = sorted(set(present) - set(expected))
    if stale:
        raise ConversionError(f"{len(stale)} family folders are not part of this conversion (e.g. "
                              f"{', '.join(stale[:5])}): a previous conversion with other settings? Not deleted; "
                              f"convert into an empty folder or remove them first")
    families.mkdir(parents=True, exist_ok=True)
    how = Counter()
    for og in expected:
        genes = groups[og]
        if tree_set == "raw_pruned":
            root = raw_tree(working, cluster_of[og], genes, internal_to_gene)
        else:
            root = resolved_tree(resolved[og], genes, prefixes)
        rows = alignment(results, alignment_source, og, cluster_of[og], genes, prefixes, gene_to_internal)
        how[write_family(paths.family(og), root, rows, release, of_version, mapping_md5)] += 1

    # 5. native/summary/: the counts and the skipped orthogroups, the split-cluster check, the model of every cluster
    split = split_check(groups, cluster_of)
    summary_dir = paths.native() / "summary"
    counts = {"orthogroups": len(groups), "families_written": len(expected),
              "orthogroups_below_min_size": len(small), "orthogroups_without_tree": len(no_tree),
              "unassigned_genes": len(unassigned), "genes": len(gene_to_internal), "species": len(species)}
    write_conversion_summary(summary_dir, counts, groups, small, no_tree, tree_set, alignment_source, min_size)
    write_split_check(summary_dir, split, groups, cluster_of, min_size)
    write_cluster_models(summary_dir, clusters, working, models_by_cluster(cache))
    say("; ".join(f"{k} {v}" for k, v in counts.items())
        + f" (families: {how['built']} built, {how['reused']} reused, {how['rebuilt']} rebuilt)")
    return {"counts": counts, "split_check": split["totals"], "tree_set": tree_set,
            "alignment_source": alignment_source, "min_family_size": min_size, "families_by_how": dict(how)}


def write_conversion_summary(folder, counts, groups, small, no_tree, tree_set, alignment_source, min_size):
    """native/summary/conversion_counts.tsv (what was converted, how many were skipped) and
    skipped_orthogroups.tsv (every orthogroup without a family folder, and why), as the family contract asks
    of a part that skips families."""
    folder.mkdir(parents=True, exist_ok=True)
    rows = [("tree_set", tree_set), ("alignment_source", alignment_source), ("min_family_size", min_size),
            *counts.items()]
    (folder / "conversion_counts.tsv").write_text("item\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in rows))
    reasons = [(og, len(groups[og]), f"fewer than {min_size} genes") for og in small] + \
              [(og, len(groups[og]), "no gene tree") for og in no_tree]
    (folder / "skipped_orthogroups.tsv").write_text(
        "orthogroup\tgenes\treason\n" + "".join(f"{og}\t{n}\t{why}\n" for og, n, why in sorted(reasons)))


def split_check(groups, cluster_of):
    """For orthogroups from a cluster split into two or more: how many, and how many species the pieces
    share (a species in two pieces of one cluster means paralogs were separated)."""
    pieces = defaultdict(list)
    for og, k in cluster_of.items():
        pieces[k].append(og)
    rows = []
    for k in sorted(pieces):
        ogs = sorted(pieces[k])
        if len(ogs) < 2:
            continue
        species_in = {og: {g.split("_", 1)[0] for g in groups[og]} for og in ogs}
        counts = Counter(s for sp in species_in.values() for s in sp)
        shared = sorted(s for s, c in counts.items() if c > 1)
        for og in ogs:
            rows.append({"cluster": f"OG{k:07d}", "orthogroup": og, "genes": len(groups[og]), "pieces": len(ogs),
                         "species": len(species_in[og]),
                         "species_shared_with_other_pieces": len(species_in[og] & set(shared))})
    totals = {"orthogroups": len(groups), "from_split_clusters": len(rows),
              "split_clusters": len({r["cluster"] for r in rows}),
              "with_a_species_in_another_piece": sum(1 for r in rows if r["species_shared_with_other_pieces"])}
    return {"rows": rows, "totals": totals}


def write_split_check(folder, split, groups, cluster_of, min_size):
    folder.mkdir(parents=True, exist_ok=True)
    cols = ["cluster", "orthogroup", "genes", "pieces", "species", "species_shared_with_other_pieces"]
    lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in split["rows"]]
    (folder / "split_clusters.tsv").write_text("\n".join(lines) + "\n")
    t = split["totals"]
    big = [r for r in split["rows"] if r["genes"] >= min_size]
    text = (f"Of {t['orthogroups']} orthogroups, {t['from_split_clusters']} come from {t['split_clusters']} "
            f"preliminary clusters that OrthoFinder split into two or more orthogroups ({len(big)} of them have "
            f"{min_size} or more genes). In {t['with_a_species_in_another_piece']} of these, at least one species "
            f"also has genes in another piece of the same cluster, so the split separated genes of one species "
            f"(paralogs) into different families; with the raw-tree-pruned tree set those genes are no longer in "
            f"one family tree. The table split_clusters.tsv lists every such orthogroup, its cluster, the number "
            f"of pieces, and how many of its species appear in another piece.\n")
    (folder / "split_clusters_summary.txt").write_text(text)


def write_cluster_models(folder, clusters, working, models):
    """cluster, sequences, model used, for every preliminary cluster with a tree. The model comes from the
    tree wrapper's record when there is one (cichlid_tree), otherwise it is 'MFP' (OrthoFinder's own IQ-TREE
    command gives no model, so ModelFinder runs) only if the cluster's IQ-TREE log says so."""
    folder.mkdir(parents=True, exist_ok=True)
    lines = ["cluster\tsequences\tmodel"]
    for k in sorted(clusters):
        cid = f"OG{k:07d}"
        if not Path(working, "Trees_ids", f"{cid}.txt").is_file():
            continue
        model = models.get(cid, (None, None))[1]
        if model is None:
            log = Path(working, "Alignments_ids", f"{cid}.log")
            text = log.read_text(errors="replace") if log.is_file() else ""
            m = re.search(r"^Command: .*? -m (\S+)", text, re.M)
            model = m.group(1) if m else ("MFP" if "ModelFinder will test" in text else "unknown")
        lines.append(f"{cid}\t{len(clusters[k])}\t{model}")
    (folder / "cluster_models.tsv").write_text("\n".join(lines) + "\n")


def run_record(run_settings, conversion, coordinator_version, catalog_md5, proteomes_md5, species, release=RELEASE):
    """The run.json record (for schemas.write_run) from the run's run_settings.json, the model setting's file
    beside it (tree_settings.json: the current rule and its history, which `set-model` changes during the run)
    and the conversion."""
    s = json.loads(Path(run_settings).read_text())
    tree_settings = Path(run_settings).with_name("tree_settings.json")
    model_file = json.loads(tree_settings.read_text()) if tree_settings.is_file() else {}
    attempts = s.get("attempts", [])
    if not attempts:
        raise ConversionError("run_settings.json has no attempts")
    first, last = attempts[0], attempts[-1]
    settings = last.get("settings", {})
    versions = last.get("versions", {})
    parameters = {
        "commands": [a.get("command") for a in attempts],
        "seed": settings.get("seed"),
        "model": settings.get("model") or settings.get("iqtree_model"),
        "model_now": model_file.get("model"), "model_since": model_file.get("since"),
        "model_history": model_file.get("history"),
        "thread_rule": settings.get("thread_rule"), "cores": settings.get("cores"),
        "tree_order": settings.get("tree_order"), "launcher_patch": settings.get("launcher_patch"),
        "tree_program": settings.get("tree_program"),
        "iqtree_argument_template": settings.get("iqtree_argument_template"),
        "versions": versions, "code": last.get("code"),
        "coordinator_versions": [{"attempt": a.get("n", i + 1), "version": a.get("coordinator_version"),
                                  **({"note": a["coordinator_version_note"]} if a.get("coordinator_version_note") else {})}
                                 for i, a in enumerate(attempts)],
        "search_reused": s.get("search_reused"),
        "species_tree": settings.get("species_tree"),
        "tree_set": conversion["tree_set"], "alignment_source": conversion["alignment_source"],
        "min_family_size": conversion["min_family_size"],
        "counts": conversion["counts"], "split_check": conversion["split_check"],
    }
    return {"release": release, "method": METHOD, "method_version": versions.get("orthofinder"),
            "parameters": parameters, "started": first.get("started"),
            "finished": last.get("finished"), "coordinator_version": coordinator_version,
            "catalog_md5": catalog_md5, "proteomes_md5": proteomes_md5, "species": list(species)}


def main(argv=None):
    p = argparse.ArgumentParser(prog="cichlid-orthofinder convert", description=__doc__.split("\n\n")[0])
    p.add_argument("--results", required=True, help="OrthoFinder's final Results_<date>/ folder")
    p.add_argument("--run-settings", required=True, help="the run's run_settings.json")
    p.add_argument("--tree-set", default="raw_pruned", choices=TREE_SETS)
    p.add_argument("--alignment", default="orthogroup", choices=ALIGNMENT_SOURCES)
    p.add_argument("--min-size", type=int, default=DEFAULT_MIN_SIZE)
    p.add_argument("--cache", help="the tree wrapper's cache (for the model of each cluster)")
    p.add_argument("--out", help="output folder (default: the coordinator's Results/<release>/OrthoFinder/)")
    p.add_argument("--release", default=RELEASE)
    args = p.parse_args(argv)
    import cichlid_data
    from cichlid_data import Catalog, get_layout, schemas
    layout = get_layout()
    catalog = Catalog.load(layout)
    species = [s.prefix for s in catalog.species(release=args.release, part=PART)]
    settings = json.loads(Path(args.run_settings).read_text())
    proteomes_md5 = catalog.proteomes_md5(species)
    if settings.get("proteomes_md5") and settings["proteomes_md5"] != proteomes_md5:
        raise SystemExit("the catalog's proteomes_md5 differs from the one the run started with; not converting")
    paths = OutputPaths(root=args.out) if args.out else OutputPaths(results_layout=layout.results(args.release, PART))
    versions = settings["attempts"][-1].get("versions", {})
    conversion = convert(args.results, paths, species, {p: layout.proteome(p).local for p in species},
                         versions.get("orthofinder") or "", args.tree_set, args.alignment, args.min_size,
                         args.release, args.cache)
    schemas.write_run(paths.run_record(), **run_record(
        args.run_settings, conversion, cichlid_data.__version__, catalog.md5(), proteomes_md5, species, args.release))
    return 0


if __name__ == "__main__":
    sys.exit(main())
