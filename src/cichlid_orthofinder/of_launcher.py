"""Start OrthoFinder with its gene trees run smallest first, without changing OrthoFinder's files.

OrthoFinder 3.1.5 has two orders for its alignment and tree jobs, "descending" (largest first, the default)
and "ascending" (smallest first); the command-line option for it is disabled in its code (process_args.py:
`self.cmd_order = "descending"  # "ascending"`, the parser branch commented out). This launcher does what
OrthoFinder's own `bin/orthofinder` does (put its `src/` on the path and call `main`), after setting that one
option on OrthoFinder's options object. The order changes only which job is handed out next; each tree is
computed from the same alignment with the same settings.

    cichlid-orthofinder-of [OrthoFinder's arguments]          (order from $CICHLID_OF_ORDER, default descending)
"""

import os
import shutil
import sys
from pathlib import Path


def orthofinder_src():
    exe = shutil.which("orthofinder")
    if not exe:
        raise SystemExit("cichlid-orthofinder-of: orthofinder is not on the PATH")
    src = Path(exe).resolve().parent / "src"
    if not (src / "orthofinder").is_dir():
        raise SystemExit(f"cichlid-orthofinder-of: OrthoFinder's source not found beside {exe}")
    return src


# The ascending order crashes in OrthoFinder 3.1.5 when a species tree is given (-s): in DoTrees the
# ascending branch computes `tree_tasksize` from `orig_commands_and_filenames`, which exists only when
# OrthoFinder infers the species tree. `tree_tasksize` is passed to RunParallelCommandsAndMoveResultsFile,
# which does not use it (only the order of the list matters), so the one-word fix cannot change any output.
# The fix is applied in memory, to the file whose MD5 is checked here; OrthoFinder's files are not changed.
TREES_MSA_MD5 = "e875a68a2dfbd6a5c390d5a28d6d8cd4"
BUG = ("            commands_and_filenames = commands_and_filenames[::-1]\n"
       "            tree_tasksize = [\n"
       "                        os.stat(item).st_size\n"
       "                        for cmd in orig_commands_and_filenames\n")
FIX = BUG.replace("for cmd in orig_commands_and_filenames", "for cmd in commands_and_filenames")


def patch_trees_msa(src):
    """Load orthofinder.tools.trees_msa with the one-line fix (in memory). -> the patched module."""
    import hashlib
    import types
    path = Path(src) / "orthofinder" / "tools" / "trees_msa.py"
    text = path.read_text()
    md5 = hashlib.md5(text.encode()).hexdigest()
    if md5 != TREES_MSA_MD5 or text.count(BUG) != 1:
        raise SystemExit(f"cichlid-orthofinder-of: {path} is not the OrthoFinder 3.1.5 file this fix was checked "
                         f"against (md5 {md5}); not running")
    import orthofinder.tools as tools_pkg
    name = "orthofinder.tools.trees_msa"
    mod = types.ModuleType(name)
    mod.__file__, mod.__package__ = str(path), "orthofinder.tools"
    sys.modules[name] = mod
    exec(compile(text.replace(BUG, FIX), str(path), "exec"), mod.__dict__)
    tools_pkg.trees_msa = mod
    return mod


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    order = os.environ.get("CICHLID_OF_ORDER", "descending")
    if order not in ("ascending", "descending"):
        raise SystemExit(f"cichlid-orthofinder-of: CICHLID_OF_ORDER must be ascending or descending, not {order}")
    src = orthofinder_src()
    sys.path.insert(0, str(src))
    patch_trees_msa(src)
    from orthofinder.run import process_args
    original = process_args.Options.__init__

    def init(self, *a, **k):
        original(self, *a, **k)
        self.cmd_order = order

    process_args.Options.__init__ = init
    from orthofinder.run.main import main as orthofinder_main
    orthofinder_main(argv)


if __name__ == "__main__":
    main()
