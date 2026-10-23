# Changes

Newest first. Entries up to `run-v1-restart` (`3dbdc2e`) name the files they changed; later ones are also
tagged commits.

## 2026-10-23: default model `full`, no approval flag (`run-v1-r2`)

- **What:** `--model` is optional and defaults to `full` (`-m MFP`, ModelFinder on every cluster) in `run`
  and `verify` (`run.DEFAULT_MODEL`, `cli.py`, `verify.py`). `hybrid` (`MFP;1000:TEST`) and other rules are
  still accepted when asked for, and `set-model` is unchanged. Every start and dry run prints the rule in
  effect and, if it is not `full`, that it departs from the sheet (`run.describe_model`). **Sheet:**
  section 4 (model).
- **What:** the approval flag is removed: `--approval`, `run.check_approval`, the requirement to give it, the
  `approval` entry in `run_settings.json` and in `run.json`'s parameters (`convert.run_record`), and its
  tests. `verify` gate 1 is now "the model setting in effect is valid" (prints the rule; PASS if it parses);
  the other gate numbers are unchanged. **Why:** the owner's decision (2026-10-23): the sheet asks that the
  owner agree to the run before it starts, not for a code gate; approvals are recorded outside the code.
- Docs: `README.md`, `docs/RUNBOOK.md`, `docs/COMPLIANCE.md`. Tests: 123 (the 3 approval tests removed; the
  2 tests that required `--model` now check the `full` default; 1 test for the new gate 1).
- **Not changed:** `treecache.py`, `of_launcher.py`, `snapshot.py`, `status.py`, and nothing in the tree
  cache key (`MFP` trees keep their keys).

## 2026-10-23: records, subcommands, species gates, runbook

- **What:** `convert.run_record` takes the model setting's history from the run's `tree_settings.json` (its
  `model`, `since` and `history`, which `set-model` changes during the run), no longer from
  `run_settings.json`, which never held it. **Sheet:** section 4 (any other model setting is recorded in
  `run.json`).
- **What:** every attempt in `run_settings.json` records `cichlid_data.__version__` (`run.py`). On the next
  real resume, an earlier attempt without it gets the version it ran under with a note: attempts started
  before 2026-10-23, such as the run's first (2026-10-22), ran under coordinator 0.13.0. A dry run writes
  nothing. `run.json` lists the version of every attempt (`parameters.coordinator_versions`). **Sheet:**
  section 5 (`run.json` records the coordinator version).
- **What:** `cichlid-orthofinder convert`, `native` and `verify` subcommands (`cli.py`), each handing its
  options to the module's own `main()`. **Why:** one command line for the whole run, for `docs/RUNBOOK.md`.
- **What:** `verify.py` gains gate 9 (every proteome re-hashed with hashlib against the catalog's per-file
  MD5), gate 10 (the part's 35 species in the catalog, in the part's list and in the coordinator's
  `config/releases.toml`, with the difference printed) and gate 11 (the MD5 of the `-s` species tree against
  the checksum the catalog recorded for it). Gates 1 and 7 say why they fail for now (no written approval
  yet; no git yet). **Sheet:** sections 3 and 9 (criterion 2).
- **What:** `.gitignore` leaves out `pilot_runs/`, `review_bundle_*/`, `*.zip`, `sample_out/`,
  `species_validation/`; `docs/RUNBOOK.md`, the command sequence from the gates to the push, with every data
  path from the coordinator.
- Tests: 125. Not changed: `treecache.py`, `of_launcher.py`, `snapshot.py`.

## 2026-10-23: the converter against the whole pack

Found by checking the code against the whole pack (`PACK-project2-orthofinder.md`: cover note,
`coordinator_spec.md` with the family contract, `species_final.md`, `releases.toml`).

- **What:** `native/summary/conversion_counts.tsv` (the settings and the counts) and
  `native/summary/skipped_orthogroups.tsv` (every orthogroup without a family folder, and why). **Why:** the
  family contract says a part that skips families says how many in its own summary under `native/summary/`;
  the counts were only in `run.json`. **Pack:** `coordinator_spec.md`, the family contract ("A family may be
  missing").
- **What:** a rerun into the same output skips every finished family whose `done.json` (inputs hash included)
  and files are unchanged, builds again a folder without `done.json`, rebuilds a family whose inputs changed,
  and refuses (without deleting) family folders the conversion would not write. Before, it refused any
  non-empty output. **Why:** the sheet and the family contract: `done.json`'s inputs hash "lets a rerun skip
  finished families", and a folder without it is unfinished and is redone. **Sheet:** sections 5 and 6.
- **What:** the output paths come from the coordinator's `ResultsLayout` (`assign(prefix)`, `family(family)`,
  `section("families")`, `native()`, `run_record()`) through `convert.OutputPaths`; a scratch folder for tests
  and samples mirrors the same structure. Before, the converter joined these paths onto the results root itself.
  **Why:** "never build a data path by hand". **Sheet:** sections 1 and 8; **pack:** `coordinator_spec.md`, Layout.
- Tests: 119 (6 more); the pilot 2 sample converted again, validated, converted twice identically, and rerun
  into the same folder (336 families reused, files unchanged).

## 2026-10-23: coordinator 0.14

- **What:** the pin moves to `cichlid-coordinator>=0.14,<0.15` (`pyproject.toml`, README), and the pin test
  expects 0.14 (`tests/test_package.py`). Coordinator 0.14.1 is installed in editable mode from a git clone
  (`~/Desktop/orthofinder_claude/cichlid-coordinator/`, `0de5afe`). **Why:** the owner's decision; the 0.13.0
  clone had been removed. The tests and `cichlid-data validate` passed on 0.14.1 first. **Sheet:** section 1
  (pin one minor version; test and validate before moving).

## 2026-10-23: the converter, the review files

- **What:** `convert.py` (OrthoFinder's output to `assign/`, `families/`, `run.json` and `native/summary/`),
  `native.py` (the curated copy of the result folder), tests `test_convert.py` and `test_native.py`, and the
  fixture `tests/fixtures/pilot2/` (five whole clusters of pilot 2's real output, 10 orthogroups, 4 unassigned
  genes, 55 proteins) with its builder `tests/fixtures/build_pilot2.py`, which replaces this machine's paths in
  IQ-TREE's logs and MCL's command line with `<local>/`.
  **Why:** the conversion was not written. **Sheet:** sections 5, 6 and 7.
- **What:** tree set, alignment source and minimum family size are parameters (`--tree-set raw_pruned|resolved`,
  `--alignment orthogroup|raw_cluster`, `--min-size`, default 4). **Why:** the owner has not decided them.
  **Sheet:** section 4 (the tree set is the owner's), section 5 (the family contract).
- **What:** `native/summary/`: `split_clusters.tsv` and `split_clusters_summary.txt` (the split-cluster check),
  `cluster_models.tsv` (cluster, sequences, model). **Why:** the check the sheet asks for before the tree set
  is decided, and a way to find the families cut from `TEST` clusters. **Sheet:** section 4.
- **What:** `docs/pilot.md`, the conversion findings on pilot 2. **Sheet:** section 9, criterion 8.

## 2026-10-23: settings, records, docs (README_for_claude_code.md, section 4)

- **What:** `--model` is required on every run (`run.py`, `cli.py`); there is no default. **Why:** the sheet's
  setting is ModelFinder on every cluster and any other needs the owner's approval; the old default was
  `hybrid`. **Sheet:** section 4 (model).
- **What:** every non-dry run (including `--resume` and `--run own_species_tree`) needs
  `--approval "<who>, <date>, <what was approved>"`, written to `run_settings.json`. **Why:** the 35-species
  run started without the owner's agreement being recorded. **Sheet:** section 4 (the owner's agreement).
- **What:** the run refuses to start if the version of OrthoFinder, IQ-TREE, FAMSA, DIAMOND or MCL cannot be read;
  FAMSA's pattern checked on this machine (`2.4.1-45c9b2b`). **Sheet:** section 4 (versions in `run.json`).
- **What:** `run_settings.json` records the package version, the SHA-256 of every source file, the tree
  program's definition, the IQ-TREE argument template, the launcher patch md5 and the `cichlid_tree` entry of
  `~/config_orthofinder_user.json` (`run.code_record()`). **Why:** the run uses `-T cichlid_tree`, not
  `-T iqtree3`; a run must say how it was made. **Sheet:** section 4. Explained in `docs/TREE_PROGRAM.md`.
- **What:** `verify.py` (not yet a subcommand): the restart gates, each printed as PASS or FAIL.
- **What:** README pin text `>=0.13,<0.14`; a note at the top of `docs/pilot.md`; the sheet in
  `docs/orthofinder_repo_spec.md` (md5 `18b33948b8cc52d7cd04e2bb56a96cfb`); `docs/SPEC_VERSION.md`;
  `docs/COMPLIANCE.md`; this file. **Sheet:** sections 1, 2 and 9.
- **Not changed:** `treecache.py`, `of_launcher.py`, `snapshot.py`, and nothing that goes into the tree cache
  key, so the stopped run's 7,619 finished trees stay valid.

## 2026-10-23: coordinator, for the review

The installed coordinator (0.13.0) ran from a clone that was removed. The owner chose to use the latest one,
0.14.1, from the GitHub zip. The review outputs (`review_bundle_2026-10-23/`) were made and validated with it,
loaded from an unpacked copy, before it was installed and the pin moved (above).
