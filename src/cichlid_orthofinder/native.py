"""What is kept under native/ (sheet section 5): a curated copy of OrthoFinder's result folder.

    cichlid-orthofinder native --results <Results_<date>/> --dest <native/our_species_tree/>

Kept: every top-level result folder and the logs; from WorkingDirectory/ only the mapping files, N0.ids.tsv,
the preliminary clusters, and the raw trees and alignments. Not copied: the search inputs and results
(Species<i>.fa, diamondDBSpecies<i>.dmnd, Blast<i>_<j>.txt.gz) and pickle/. They stay, read-only, on the
machine that ran them, reusable with -b.

A copied OrthoFinder folder still points at the original through absolute paths in its Log.txt, so the copy
is a record, not a folder OrthoFinder can be restarted from.
"""

import argparse
import shutil
import sys
from pathlib import Path

WORKING_KEEP = ("SpeciesIDs.txt", "SequenceIDs.txt", "N0.ids.tsv", "clusters_OrthoFinder_I1.2.txt", "Trees_ids",
                "Alignments_ids")
NOT_COPIED = ("Species*.fa", "diamondDBSpecies*.dmnd", "Blast*_*.txt.gz", "pickle")


def keep_list(results):
    """-> [(path relative to the result folder, is it a folder)] of what the curated copy holds, sorted."""
    results = Path(results)
    keep = []
    for entry in sorted(results.iterdir()):
        if entry.name == "WorkingDirectory":
            for name in WORKING_KEEP:
                if (entry / name).exists():
                    keep.append((f"WorkingDirectory/{name}", (entry / name).is_dir()))
        else:
            keep.append((entry.name, entry.is_dir()))
    return keep


def size(path):
    path = Path(path)
    if path.is_file():
        return path.stat().st_size
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def curate(results, dest):
    """Copy the keep-list of `results` into `dest` (a new folder). -> {top-level folder or file: bytes}."""
    results, dest = Path(results), Path(dest)
    if dest.exists():
        raise FileExistsError(f"{dest} exists: the curated copy goes into a new folder")
    for rel, is_dir in keep_list(results):
        src, target = results / rel, dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if is_dir:
            shutil.copytree(src, target)
        else:
            shutil.copy2(src, target)
    return {e.name: size(e) for e in sorted(dest.iterdir())}


def main(argv=None):
    p = argparse.ArgumentParser(prog="cichlid-orthofinder native", description=__doc__.split("\n\n")[0])
    p.add_argument("--results", required=True, help="OrthoFinder's Results_<date>/ folder")
    p.add_argument("--dest", required=True, help="the new curated folder, e.g. native/our_species_tree/")
    args = p.parse_args(argv)
    sizes = curate(args.results, args.dest)
    for name, n in sizes.items():
        print(f"{n / 1e6:12.1f} MB  {name}")
    print(f"{sum(sizes.values()) / 1e6:12.1f} MB  total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
