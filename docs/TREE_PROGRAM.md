# The tree program: `-T cichlid_tree` instead of `-T iqtree3`

The sheet (section 4) asks for `-T iqtree3`. The run uses `-T cichlid_tree`: a tree program registered in
OrthoFinder's user configuration, `~/config_orthofinder_user.json`, by `run.register_tree_program()`:

```json
"cichlid_tree": {
  "program_type": "tree",
  "cmd_line": "cichlid-orthofinder-tree INPUT PATH/IDENTIFIER",
  "output_filename": "PATH/IDENTIFIER.treefile"
}
```

`cichlid-orthofinder-tree` (`treecache.py`) runs the same `iqtree3` on the same alignment, with OrthoFinder's
`-s` and `--prefix` and IQ-TREE's `-quiet`, and adds three things:

```
iqtree3 -s <alignment> --prefix <prefix> -m <model> -seed <seed> -T <threads> -quiet
```

- `-seed <seed>`: the fixed seed the sheet requires. OrthoFinder's built-in command has no seed, and in the
  pilot a rerun changed 108 of 129 trees. This is why the wrapper exists.
- `-m <model>`: `MFP` (ModelFinder, what the built-in command does by default), written out, unless the
  owner approves another setting by size (the `--model` rule; `hybrid` = `MFP;1000:TEST`).
- `-T <threads>`: threads by cluster size within a shared core budget (`--thread-rule`), which is the
  machine's business, not science.

It also keeps a cache keyed by a SHA-256 of these settings and the alignment, and IQ-TREE's checkpoints, so
a stopped run resumes without rebuilding finished trees.

`run_settings.json` records the tree program's definition (`settings.tree_program`), the IQ-TREE argument
template (`settings.iqtree_argument_template`), the `cichlid_tree` entry actually found in
`~/config_orthofinder_user.json` (`code.orthofinder_user_config_entry`), the package version, and the
SHA-256 of every source file. `run.json` carries them into its `parameters`.

The same configuration file also holds two entries from the pilot's tests, `iqtree3_mt` and
`iqtree3_mt_seed`. The run does not use them.
