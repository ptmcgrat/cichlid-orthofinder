"""Interim results while the largest trees are still running: an OrthoFinder output covering only the
clusters whose trees are finished, made from a copy so the running run is untouched.

OrthoFinder's last steps (`-fgt`, "from gene trees") take the trees in the working folder and do the rest:
the split into final orthogroups, the orthologues, the duplications and the result files. Clusters without
a tree are left out (their genes are reported unassigned); every other orthogroup comes out with the same
genes and the same resolved tree as in the final run (tested, docs/pilot.md). OrthoFinder numbers the
orthogroups by size, so the IDs differ from the final run's: `mapping()` writes the correspondence.

Three things the copy must get right:
- `-fgt` finds its files through absolute paths written in Log.txt (WorkingDirectory_Base,
  WorkingDirectory_Trees, FN_Orthogroups). A copy keeps pointing at the original unless they are rewritten,
  and then silently uses the original's trees. The copy's Log.txt gets all three, pointing at the copy.
- OrthoFinder writes WorkingDirectory_Trees to Log.txt only after the tree stage, so a running run's log does
  not have it: the trees are read from the current Results folder's WorkingDirectory instead.
- after a resume, the species files are in the first working folder (WorkingDirectory_Base) and the
  clusters, alignments and trees in the current one; both are read.
"""

import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

LOG_KEYS = ("WorkingDirectory_Base", "WorkingDirectory_Trees", "FN_Orthogroups")
SKIP = re.compile(r"^(Blast\d+_\d+\.txt(\.gz)?|diamondDBSpecies\d+\.dmnd)$")   # the search: not needed


def log_paths(results):
    """-> {key: path} for the paths Log.txt records so far (the last value of each)."""
    out = {}
    for line in Path(results, "Log.txt").read_text(errors="replace").splitlines():
        for k in LOG_KEYS:
            if line.startswith(k + ":"):
                out[k] = line.split(":", 1)[1].strip()
    return out


def _copy_working(base, current, dest):
    """Copy what -fgt needs into dest: the species files from the first working folder (not the search),
    then the current working folder's files (clusters, sequences, pickles), its finished trees and their
    alignments. -> number of trees copied."""
    dest.mkdir(parents=True)
    for src in dict.fromkeys([Path(base), Path(current)]):          # base first; the current one wins
        for p in src.iterdir():
            if p.is_file() and not SKIP.match(p.name):
                shutil.copy2(p, dest / p.name)
            elif p.is_dir() and p.name in ("Sequences_ids", "dependencies", "pickle"):
                if (dest / p.name).exists():
                    shutil.rmtree(dest / p.name)
                shutil.copytree(p, dest / p.name)
    (dest / "Trees_ids").mkdir()
    (dest / "Alignments_ids").mkdir()
    n = 0
    for t in sorted(Path(current, "Trees_ids").glob("*.txt")):
        shutil.copy2(t, dest / "Trees_ids" / t.name)       # OrthoFinder puts a tree here only when it is complete
        aln = Path(current, "Alignments_ids", t.stem + ".fa")
        if aln.is_file():
            shutil.copy2(aln, dest / "Alignments_ids" / aln.name)
        n += 1
    return n


def _rewrite_log(src, dst, wd):
    """Copy Log.txt with the three paths pointing at the copy's working folder (added if missing)."""
    wd = Path(wd)
    clusters = sorted(wd.glob("clusters_OrthoFinder_I*.txt_id_pairs.txt"))
    if not clusters:
        raise RuntimeError(f"no clusters_OrthoFinder_I*.txt_id_pairs.txt in {wd}")
    keep = [l for l in Path(src, "Log.txt").read_text(errors="replace").splitlines()
            if not any(l.startswith(k + ":") for k in LOG_KEYS)]
    keep += [f"WorkingDirectory_Base: {wd}/", f"FN_Orthogroups: {clusters[-1]}", f"WorkingDirectory_Trees: {wd}/"]
    Path(dst, "Log.txt").write_text("\n".join(keep) + "\n")


def _count(p):
    d = Path(p).read_bytes()
    return d.count(b"\n>") + d.startswith(b">")


def take(results, dest, species_tree=None, threads=8, runner=subprocess.run, say=print):
    """Snapshot the run whose current Results_* folder is `results` into `dest` (a new folder).
    -> (the snapshot's results folder, its manifest)."""
    results, dest = Path(results).resolve(), Path(dest).resolve()
    if dest.exists():
        raise RuntimeError(f"snapshot folder exists: {dest}")
    current = results / "WorkingDirectory"
    base = Path(log_paths(results).get("WorkingDirectory_Base") or current)
    stage = dest / "Results_snapshot"
    wd = stage / "WorkingDirectory"
    n_trees = _copy_working(base, current, wd)
    _rewrite_log(results, stage, wd)

    have = {p.stem for p in (wd / "Trees_ids").glob("*.txt")}
    excluded = sorted(((p.stem, _count(p)) for p in (wd / "Sequences_ids").glob("*.fa")
                       if p.stem not in have and _count(p) >= 3), key=lambda x: -x[1])
    cmd = ["orthofinder", "-fgt", str(stage), "-M", "msa", "-T", "cichlid_tree", "-t", str(threads)]
    if species_tree:
        cmd += ["-s", str(species_tree)]
    say("$ " + " ".join(cmd))
    with open(dest / "orthofinder.log", "w") as log:
        rc = runner(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    # -fgt writes a new Results_<date> folder beside the staged one
    made = sorted((p for p in dest.glob("Results_*") if p != stage and (p / "Orthogroups").is_dir()),
                  key=lambda p: p.stat().st_mtime)
    if rc != 0 or not made:
        raise RuntimeError(f"orthofinder -fgt failed (exit {rc}); see {dest / 'orthofinder.log'}")
    manifest = {
        "taken": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source_results": str(results), "trees_included": n_trees,
        "clusters_excluded": len(excluded), "genes_excluded": sum(n for _, n in excluded),
        "excluded": [{"cluster": c, "sequences": n} for c, n in excluded],
        "snapshot_results": str(made[-1]), "command": cmd,
        "note": "orthogroup IDs are OrthoFinder's numbering for this snapshot and differ from the final run's; "
                "see mapping.tsv once the final run is done",
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    say(f"snapshot: {n_trees} trees, {len(excluded)} clusters left out ({manifest['genes_excluded']} genes); "
        f"results in {made[-1]}")
    return made[-1], manifest


def ready(results, size):
    """True when every cluster of 3 to `size` sequences has its tree and some larger cluster has not: the
    moment for interim results up to that size. False while the tree stage has not begun."""
    wd = Path(results) / "WorkingDirectory"
    seqs, trees = wd / "Sequences_ids", wd / "Trees_ids"
    if not seqs.is_dir() or not trees.is_dir():
        return False
    have = {p.stem for p in trees.glob("*.txt")}
    large_left = 0
    for p in seqs.glob("*.fa"):
        if p.stem in have:
            continue
        n = _count(p)
        if 3 <= n <= size:
            return False                                   # a cluster of this size or smaller has no tree yet
        if n > size:
            large_left += 1
    return large_left > 0


class Watcher:
    """Interim results at several sizes: once every cluster of up to `size` sequences has its tree (for each
    size in `sizes`, smallest first), a snapshot; when the run has finished, a mapping from each snapshot to
    the final results. Call poll() regularly while OrthoFinder runs and finish() when it has completed.
    `state` (kept in run_settings.json by the caller) holds the snapshots taken, so a resumed run does not
    take them again. If a larger size is already ready when a smaller one was missed (a resume, say), only
    the larger is taken."""

    def __init__(self, work, output_dir, sizes, species_tree, threads=8, state=None, say=print):
        self.work, self.output_dir = Path(work), Path(output_dir)
        self.sizes = sorted({int(s) for s in ([sizes] if isinstance(sizes, int) else sizes)})
        self.species_tree, self.threads, self.say = species_tree, threads, say
        self.state = state if state is not None else {}
        self.state.setdefault("snapshots", [])

    def taken(self, size):
        return any(s["size"] == size for s in self.state["snapshots"])

    def poll(self):
        from .run import latest_results
        pending = [s for s in self.sizes if not self.taken(s) and s > max([x["size"] for x in self.state["snapshots"]],
                                                                            default=0)]
        if not pending:
            return False
        try:
            results = latest_results(self.output_dir)
        except (RuntimeError, FileNotFoundError):
            return False
        due = [s for s in pending if ready(results, s)]
        if not due:
            return False
        size = max(due)                                   # skip smaller sizes overtaken while not watching
        dest = self.work / f"snapshot_upto{size}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.say(f"interim results: every cluster of up to {size} sequences has its tree; snapshot in {dest}")
        try:
            snap, manifest = take(results, dest, self.species_tree, self.threads, say=self.say)
        except Exception as e:                            # recorded; the run goes on; retried at the next poll
            self.state.setdefault("errors", []).append(f"{datetime.now().isoformat(timespec='seconds')} "
                                                       f"up to {size}: {e}")
            self.say(f"interim snapshot up to {size} failed (the run continues; retried later): {e}")
            if dest.exists():
                shutil.rmtree(dest, ignore_errors=True)
            return True                                   # state changed: the caller saves it
        self.state["snapshots"].append(dict(size=size, folder=str(dest), results=str(snap), taken=manifest["taken"],
                                            clusters_excluded=manifest["clusters_excluded"],
                                            genes_excluded=manifest["genes_excluded"],
                                            skipped_smaller=[s for s in due if s < size]))
        return True

    def finish(self, final_results):
        out = []
        for s in self.state["snapshots"]:
            path = self.work / f"snapshot_upto{s['size']}_mapping.tsv"
            summary = mapping(s["results"], final_results, path)
            s.update(mapping=str(path), **{f"mapping_{k}": v for k, v in summary.items()})
            self.say(f"snapshot up to {s['size']}: {summary['matched']} of {summary['orthogroups']} interim "
                     f"orthogroups found in the final results, {summary['same_tree']} with the same tree; {path}")
            out.append(summary)
        return out


def _orthogroups(results):
    out = {}
    for line in Path(results, "Orthogroups", "Orthogroups.txt").read_text().splitlines():
        og, _, genes = line.partition(":")
        out[og] = frozenset(genes.split())
    return out


def _resolved(results):
    f = Path(results, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    return dict(l.split(": ", 1) for l in f.read_text().splitlines() if ": " in l)


def mapping(snapshot_results, final_results, out):
    """Write snapshot orthogroup -> final orthogroup (matched by identical genes), and whether the resolved
    tree is identical. Snapshot orthogroups of one gene (genes left out) are skipped. -> summary dict."""
    S, F = _orthogroups(snapshot_results), _orthogroups(final_results)
    by_genes = {genes: og for og, genes in F.items()}
    RS, RF = _resolved(snapshot_results), _resolved(final_results)
    rows, matched, same_tree = [], 0, 0
    for og, genes in sorted(S.items()):
        if len(genes) < 2:
            continue
        final = by_genes.get(genes, "")
        tree = bool(final) and og in RS and RS.get(og) == RF.get(final)
        matched += bool(final)
        same_tree += tree
        rows.append(f"{og}\t{final}\t{len(genes)}\t{'yes' if tree else ('' if og not in RS else 'no')}\n")
    Path(out).write_text("snapshot_og\tfinal_og\tgenes\tsame_resolved_tree\n" + "".join(rows))
    return {"orthogroups": len(rows), "matched": matched, "same_tree": same_tree}
