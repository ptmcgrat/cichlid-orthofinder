#!/usr/bin/env python3
"""A stand-in for iqtree3, for offline tests of the tree wrapper.

Behaves like IQ-TREE where the wrapper cares: `--version`; with `-s ALN --prefix P ... -seed S` it writes
P.log, P.ckp.gz (its checkpoint), then P.treefile and P.iqtree, and FINISHES the log. A rerun on the same
prefix with a checkpoint "resumes" (logs CHECKPOINT). Controlled by environment variables:

    FAKE_IQTREE_CALLS   file: one line appended per tree run ("<prefix> fresh|resume <model>")
    FAKE_IQTREE_CRASH   "checkpoint": exit after writing the checkpoint (a kill during the run)
                        "after_tree": exit after finishing, as if killed before the wrapper recorded it
                        "fail": exit 2 at once
                        "bad_resume": fail if a checkpoint exists (an unusable checkpoint)
The tree depends only on the alignment's names and the seed, so identical inputs give identical trees.
"""
import hashlib
import os
import sys
from pathlib import Path

args = sys.argv[1:]
if args == ["--version"]:
    print("IQ-TREE version 9.9.9 for testing")
    sys.exit(0)
aln = Path(args[args.index("-s") + 1])
prefix = args[args.index("--prefix") + 1]
seed = args[args.index("-seed") + 1]
crash = os.environ.get("FAKE_IQTREE_CRASH", "")
ckp = Path(prefix + ".ckp.gz")
resuming = ckp.exists()
if os.environ.get("FAKE_IQTREE_CALLS"):
    with open(os.environ["FAKE_IQTREE_CALLS"], "a") as fh:
        model = args[args.index("-m") + 1] if "-m" in args else "-"
        fh.write(f"{prefix} {'resume' if resuming else 'fresh'} {model}\n")
if crash == "fail" or (crash == "bad_resume" and resuming):
    sys.exit(2)
with open(prefix + ".log", "a") as log:
    if resuming:
        log.write(f"CHECKPOINT: Resuming analysis from {ckp}\n")
    log.write("ModelFinder will test up to 1232 protein models\n")
    ckp.write_text("checkpoint\n")
    if crash == "checkpoint":
        sys.exit(137)
    names = sorted(l[1:].split()[0] for l in aln.read_text().splitlines() if l.startswith(">"))
    h = hashlib.sha256(("|".join(names) + seed).encode()).hexdigest()[:8]
    Path(prefix + ".treefile").write_text("(" + ",".join(names) + f"):{int(h, 16) % 7};\n")
    Path(prefix + ".iqtree").write_text("report\n")
    log.write("Total wall-clock time used: 0.01 sec\n")
if crash == "after_tree":
    sys.exit(137)
