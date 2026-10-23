"""Build tests/fixtures/pilot2/ from pilot 2's real OrthoFinder output (5 species, 2,183 proteins).

    python tests/fixtures/build_pilot2.py <pilot 2 Results_<date>/> <pilot 2 input folder with <prefix>.faa> <dest>

Whole preliminary clusters are kept (so every raw tree and alignment is the real one), with every orthogroup
cut from them, a few unassigned genes, and the matching lines of every table and mapping file. The proteomes
hold exactly those genes, in their original order, so the mapping checks of the converter hold.
"""
import re
import shutil
import sys
from pathlib import Path

# this machine's paths in IQ-TREE's logs and MCL's command line are replaced by a placeholder
LOCAL_PATH = re.compile(r"(?:/tmp/claude-[^\s\"']*?/scratchpad|/home/[^/\s\"']+)/")


def scrub(path):
    path.write_text(LOCAL_PATH.sub("<local>/", path.read_text(errors="replace")))

CLUSTERS = (0, 5, 52, 79, 316)        # split in 3; split in 2; three small orthogroups; one whole; the None alignment
UNASSIGNED = 3                         # unassigned genes kept (the first of each of the first species)


def main(results, inputs, dest):
    results, inputs, dest = Path(results), Path(inputs), Path(dest)
    W, out = results / "WorkingDirectory", dest / "Results_Oct20"
    if dest.exists():
        shutil.rmtree(dest)
    (out / "WorkingDirectory").mkdir(parents=True)

    pairs = (W / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").read_text()
    head, body = pairs.split("begin", 1)
    body = body.rsplit(")", 1)[0]
    rows = {int(c.split()[0]): c.split()[1:] for c in body.split("$") if c.split()}
    keep_internal = {i for k in CLUSTERS for i in rows[k]}
    seqids = {l.split(":")[0]: l for l in (W / "SequenceIDs.txt").read_text().splitlines() if l.strip()}
    gene_of = {k: v.split(":", 1)[1].strip() for k, v in seqids.items()}

    unassigned_lines, kept_unassigned, extra = [], set(), 0
    in_clusters = {gene_of[i] for i in keep_internal}
    lines = (results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").read_text().splitlines()
    for line in lines[1:]:
        cells = line.split("\t")
        genes = {c for c in cells[1:] if c}
        if genes & in_clusters or (extra < UNASSIGNED and genes):   # every one in a kept cluster, and a few more
            extra += not (genes & in_clusters)
            unassigned_lines.append(line)
            kept_unassigned |= genes
    keep_genes = {gene_of[i] for i in keep_internal} | kept_unassigned
    keep_internal |= {i for i, g in gene_of.items() if g in kept_unassigned}

    # mapping files and clusters
    shutil.copy2(W / "SpeciesIDs.txt", out / "WorkingDirectory" / "SpeciesIDs.txt")
    (out / "WorkingDirectory" / "SequenceIDs.txt").write_text(
        "".join(seqids[i] + "\n" for i in seqids if i in keep_internal))
    text = "".join(f"{k}      {' '.join(rows[k])} $\n" for k in sorted(CLUSTERS))
    (out / "WorkingDirectory" / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text(head + "begin\n" + text + ")\n")
    shutil.copy2(W / "clusters_OrthoFinder_I1.2.txt", out / "WorkingDirectory" / "clusters_OrthoFinder_I1.2.txt")
    for folder, suffixes in (("Trees_ids", (".txt",)), ("Alignments_ids", (".fa", ".log"))):
        (out / "WorkingDirectory" / folder).mkdir()
        for k in CLUSTERS:
            for s in suffixes:
                f = W / folder / f"OG{k:07d}{s}"
                if f.exists():
                    shutil.copy2(f, out / "WorkingDirectory" / folder / f.name)

    # orthogroups, alignments, resolved trees
    (out / "Orthogroups").mkdir()
    ogs = []
    for line in (results / "Orthogroups" / "Orthogroups.txt").read_text().splitlines():
        og, genes = line.split(":", 1)
        if set(genes.split()) <= keep_genes:
            ogs.append(line)
    (out / "Orthogroups" / "Orthogroups.txt").write_text("".join(l + "\n" for l in ogs))
    (out / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(lines[0] + "\n" + "".join(l + "\n" for l in unassigned_lines))
    names = {l.split(":")[0] for l in ogs}
    (out / "MultipleSequenceAlignments").mkdir()
    for og in sorted(names):
        f = results / "MultipleSequenceAlignments" / f"{og}.fa"
        if f.exists():
            shutil.copy2(f, out / "MultipleSequenceAlignments" / f.name)
    (out / "Resolved_Gene_Trees").mkdir()
    resolved = [l for l in (results / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").read_text().splitlines()
                if l.split(":")[0] in names]
    (out / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text("".join(l + "\n" for l in resolved))
    log = (results / "Log.txt").read_text().splitlines()
    (out / "Log.txt").write_text(log[0] + "\n")

    # the proteomes: exactly the kept genes, in their original order
    (dest / "proteomes").mkdir()
    for faa in sorted(inputs.glob("*.faa")):
        keep, out_lines = False, []
        for line in faa.read_text().splitlines():
            if line.startswith(">"):
                keep = line[1:].split()[0] in keep_genes
            if keep:
                out_lines.append(line)
        (dest / "proteomes" / faa.name).write_text("".join(l + "\n" for l in out_lines))
    for f in [*out.rglob("*.log"), *out.glob("WorkingDirectory/clusters_*")]:
        scrub(f)
    print(f"{len(names)} orthogroups, {len(keep_genes)} genes, {len(CLUSTERS)} clusters -> {dest}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
