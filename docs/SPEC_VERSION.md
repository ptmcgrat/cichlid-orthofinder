# Spec version

The instruction sheet and coordinator this repository is written against.

| | |
|---|---|
| Instruction sheet | `docs/orthofinder_repo_spec.md`, from the pack `PACK-project2-orthofinder.md` |
| Dispatch | Dispatch 5, 2026-10-22 |
| Sheet md5 | `18b33948b8cc52d7cd04e2bb56a96cfb` (extracted from the pack on 2026-10-23; identical to the Dispatch 5 original) |
| Coordinator version | 0.14.1 (pin `cichlid-coordinator>=0.14,<0.15` in `pyproject.toml`). The sheet and the pack were written against 0.13.0 |
| Move to 0.14 | 2026-10-23, the owner's decision: the 0.13.0 clone had been removed. 0.14.1 (GitHub `main`, `0de5afe`) is installed in editable mode from a git clone at `~/Desktop/orthofinder_claude/cichlid-coordinator/`. The tests and `cichlid-data validate` passed on 0.14.1 first, as the sheet (section 1) requires before moving the pin. The coordinator spec's changes from 0.13.0 to 0.14.1 (reference names for the homology project, communication rules) touch no OrthoFinder format; the species list and the species tree are unchanged. |
| Date | 2026-10-23 |
