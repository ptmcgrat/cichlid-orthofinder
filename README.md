# cichlid-orthofinder

This repository is project 2 of the cichlid protein phylogeny resource. It runs OrthoFinder on the 35 proteomes of
release v1 and turns OrthoFinder's results into the coordinator's format: assignment files, family folders and a
run record, under `Results/v1/OrthoFinder/`. The next project (the homology decider) reads them from there.

It follows the coordinator's instruction sheet for this repository, `docs/orthofinder_repo_spec.md` in the
coordinator (a verbatim copy is in `docs/` here), and works with coordinator version 0.14.

**Status (2026-10-24):**
- **The code is finished and tested** (tag `run-v1-r6`: early results now also come out in the coordinator's format,
  and results of a resumed run convert correctly; see `CHANGES.md`).
- **The input species and the output format match the coordinator,** and every check passes.
- **The full 35-species run is paused.** It waits for the coordinator to approve the run-time choices listed in
  "Waiting for coordinator approval" below. Then it restarts, reusing the 7,619 gene trees already finished.

---

## What it does, step by step

1. **Gets the inputs from the coordinator.**
   - It runs `cichlid-data check --release v1 --part OrthoFinder --fetch` and stops if that fails.
   - It reads the list of 35 species, their proteome files and the species tree from the coordinator's catalog.
     Nothing is typed by hand, and no copy of the species list or the tree is kept in this repository.
   - OrthoFinder reads links to the proteomes from a separate scratch folder, so nothing is ever written into the
     coordinator's input folders.
2. **Runs OrthoFinder 3.1.5** with the settings the sheet agrees (section 4): `-M msa`, IQ-TREE for the gene trees,
   and our species tree with `-s`. OrthoFinder:
   - compares every protein with every other (DIAMOND);
   - groups them into clusters (MCL);
   - aligns each cluster (FAMSA) and builds its gene tree (IQ-TREE 3.1.3, with ModelFinder choosing the model for
     every tree, as the sheet asks);
   - uses the species tree to split the clusters into orthogroups.
3. **Converts OrthoFinder's results** into the coordinator's format (details in the next section):
   - one assignment file per species;
   - one family folder per orthogroup;
   - a run record;
   - a trimmed copy of OrthoFinder's own result folder.
4. **Checks the result** with the coordinator's own validators (`cichlid-data validate`) before anything is handed over.

Every setting, program version and command is recorded, so the run can always say exactly how it was made.

---

## What it writes, and how that matches the coordinator

All output goes to the folder the coordinator gives for this part (`cichlid-data where results v1 OrthoFinder`). The
structure is the coordinator's (`coordinator_spec.md`, the Layout section, and sheet section 5):

```
Results/v1/OrthoFinder/
├── assign/<prefix>.assign.tsv     one per species (Mzeb.assign.tsv, Acal.assign.tsv, ...)
├── families/<orthogroup>/         one folder per orthogroup (OG0000000/, OG0000001/, ...)
│   ├── tree.nwk                   the gene tree
│   ├── combined.aln               the alignment
│   └── done.json                  written last; marks the family as finished
├── native/
│   ├── our_species_tree/          a trimmed copy of OrthoFinder's own result folder
│   └── summary/                   counts, skipped orthogroups and why, the split-cluster check
└── run.json                       the run record
```

`names/`, `orthology/` and `residues/` are left empty on purpose. The coordinator gives them to the homology decider
(project 3) and the analysis project (project 4). The sheet says this project makes no orthology calls and no names
(section 8).

| Output | What it holds | Where it comes from | How it matches the coordinator |
|---|---|---|---|
| `assign/<prefix>.assign.tsv` | Every gene of the species, in proteome order: its gene ID, protein length, and orthogroup (empty for a gene in none) | The gene IDs and lengths come from the coordinator's proteome; the orthogroup from OrthoFinder's `Orthogroups.txt` | The 16 columns in the coordinator's order, written with the coordinator's own writer (`schemas.ASSIGN.write`). The score, coverage and fusion columns stay empty, as the sheet says |
| `families/<OG>/tree.nwk` | The orthogroup's gene tree | OrthoFinder's tree for the cluster, cut down to the orthogroup's genes | Every leaf is renamed back to our gene ID (`Mzeb_101487032`). Branch lengths are kept, there are no support values, and the tree is unrooted with `"rooted": false`, as the sheet says |
| `families/<OG>/combined.aln` | The orthogroup's alignment | OrthoFinder's alignment | Rows renamed to the same gene IDs. A family is refused unless there is exactly one row for every leaf (pilot 2 found broken alignments) |
| `families/<OG>/done.json` | `release`, `rooted`, `tree_mode` (`denovo`) and an inputs hash | Written by the conversion, last | Exactly the keys of the family contract. A rerun skips families whose inputs did not change |
| `native/our_species_tree/` | OrthoFinder's result folders and logs | Copied from OrthoFinder's output | Exactly the sheet's keep-list (section 5), including the six files and folders it names from `WorkingDirectory/`. The 8.7 GB of search files are left out. Folder sizes are printed before any push |
| `run.json` | Release, method and version, every setting (seed, model, threads), the commands, the coordinator version, the catalog and proteome checksums, the 35 species | This run's own record | All required keys, written with the coordinator's own writer (`schemas.write_run`) |

Orthogroups of fewer than 4 genes get no family folder (the minimum size can be changed). Their genes still appear in
`assign/`, and they are listed in `native/summary/skipped_orthogroups.tsv`, as the family contract allows.

---

## Compliance with the coordinator

| Coordinator's requirement | Status | How it was checked |
|---|---|---|
| Exactly the 35 species of v1, in the catalog's order | Met | The coordinator's config, the catalog and the code give the same 35 in the same order |
| The coordinator's proteomes, unchanged | Met | Every file's checksum matches the catalog's. The combined checksum is `5df42f6efb77659a03c6a82d6d9c7b26` and there are 848,665 proteins, the coordinator's values |
| The catalog's species tree, byte for byte | Met | The tree's checksum, `b1296181e36810a6607fcf2988f786c6`, is the same in the coordinator's config, the catalog and the run |
| `cichlid-data check` passes before a run | Met | 36 of 36 inputs ready. The run refuses to start otherwise |
| Every output passes `cichlid-data validate` | Met on test outputs | Assign files, families and `run.json` all validate. The full run's output is validated the same way before hand-over |
| Never type a data path; no rclone; never write into an input folder | Met | A code search finds none; all paths come from the coordinator's layout. The input folders were listed before and after the test runs, and nothing changed |
| A fixed `-seed`, recorded in `run.json` | Met | Seed 1, recorded |
| ModelFinder on every cluster | Met | It is the default (`-m MFP`); every start prints the model in effect |
| The program versions of the sheet (OrthoFinder 3.1.5, IQ-TREE 3.1.3, FAMSA 2.4.1, DIAMOND 2.2.8, MCL 22-282) | Met | Read from each program at every start and recorded. The run refuses to start if one can't be read |
| Pinned to one tested coordinator version | Met | `cichlid-coordinator>=0.14,<0.15`, tested on 0.14.1. The coordinator's rule allows moving the pin after testing. The coordinator's handoff table still lists this project on 0.13 |

---

## How we know the species are right

1. **One source.** The species list, the protein files and the species tree come only from the coordinator's
   catalog, at every run. This repository holds no copy of them.
2. **The list is checked three ways.** The coordinator's config (`releases.toml`), its catalog, and what the code
   reads give the same 35 species in the same order.
3. **Every file is checked by checksum.** Each protein file's checksum matches the one the catalog recorded. The
   combined checksum of all 35 (`5df42f6e…`) and the protein count (848,665) match the coordinator's values.
4. **The tree is checked the same way.** It is the coordinator's file byte for byte (`b1296181…`).
5. **The run refuses to start if anything differs.** Tested: changing one letter in a copied protein file, or using
   a different tree, stops it.
6. **The saved comparison was checked before reuse.** The protein-against-protein comparison from earlier was
   compared, sequence by sequence, with the 35 protein files. It is reused only because it matches them exactly.

---

## What was tested

Everything below passed on 2026-10-23 and was rechecked on 2026-10-24.

**The basics:**
- **129 offline tests** run in about 2 seconds: `PYTHONPATH=src python -m unittest discover tests`.
- **Restart gates, 11 of 11:** `cichlid-orthofinder verify` checks the inputs, the tree, the versions, the data folder
  and the code before a restart.

**Small real runs** (5 species, 306 proteins cut from the real proteomes; scratch folders only):
- **End to end:** run, convert, then the coordinator's validators. Everything validates.
- **Conversion:** all four combinations of tree set and alignment source validate, and converting twice gives
  identical files.
- **Crash recovery:**
  - A run was killed at different points: during the search, while choosing a model, during a tree search, and
    after some trees. Each time, the restart reused every finished tree and continued the unfinished ones.
  - A run was also finished by hand after a crash.
  - Single-threaded results were identical to an uninterrupted run.
- **Missing trees are caught:** when tree jobs fail, OrthoFinder still reports success. The run notices the missing
  trees and refuses to accept the result until they are rebuilt.
- **Core limit:** the total number of threads in use never went above the limit given (13 runs, measured live).
- **Job order:** running the smallest trees first leaves OrthoFinder's output byte-identical (139 files, 33 trees).

**The real data:**
- **Species:** the checks in "Compliance" above, run on the real 35 proteomes and the real species tree.


---

## Features (how the run is made robust)

None of these change the inputs or the output format.

- **Restart without losing work.**
  - Each gene tree goes through a small program of ours (`cichlid-orthofinder-tree`). It runs the same IQ-TREE
    command OrthoFinder would, adds the fixed seed, and keeps every finished tree.
  - After a crash, `--resume` reuses all finished trees and continues the unfinished ones from IQ-TREE's own
    checkpoint files.
  - A search cut short by a crash is detected and run again.
- **A core limit.** Large trees get more threads (default: 1 thread per tree, 2 from 600 sequences, 8 from 1,500,
  16 from 3,000). All trees together never use more than the cores given.
- **Smallest trees first,** so early results arrive sooner. This uses an option already inside OrthoFinder, plus a
  one-word fix applied in memory; OrthoFinder's files are not changed.
- **Early results.** Interim snapshots are taken once every family up to 300 and up to 1,000 genes has its tree.
  Each one is converted straight away into the coordinator's format (same `assign/`, `families/`, `native/summary/`
  and `run.json`, same writers and checks), inside the snapshot's own folder:
  `<run folder>/snapshot_upto<N>_<time>/Results/v1/OrthoFinder/`. It never goes into the coordinator's Results
  folder and is not pushed. Two differences from the final results, both written in its `run.json`
  (`parameters.interim`): the orthogroup numbers are the snapshot's own (a `mapping.tsv` links them to the final ones
  when the run ends), and genes of families still running have no family yet. Tested on a small real run: all 31
  interim families were byte for byte the same as their final ones.
- **Progress report.** `cichlid-orthofinder status <run folder>` shows the stage, the trees done and the expected dates.
- **The model by family size** can be set: the default is ModelFinder for all, and changing it needs approval.
  It can even be changed while the run is going (`set-model`).

---

## Waiting for coordinator approval

The inputs and the output format are fixed by the coordinator and are not open. These are choices about **how the
trees are made**.

### The exact commands we run

**OrthoFinder** (started by `cichlid-orthofinder run`; printed by `--dry-run` and recorded in `run_settings.json`):

```
cichlid-orthofinder-of -f <input links> -o <output> -M msa -T cichlid_tree -s <coordinator's species tree> -t 32
            (on a restart: -b <the finished search> instead of -f/-o)
```

`cichlid-orthofinder-of` is OrthoFinder 3.1.5 itself, started with its smallest-first job order switched on
(see "Job order" below); the options are OrthoFinder's own.

| Flag | Meaning | The sheet (section 4) | Ours |
|---|---|---|---|
| `-M msa` | Alignments and maximum-likelihood gene trees | `-M msa` | same |
| `-T iqtree3` / `-T cichlid_tree` | Which program builds the trees | `-T iqtree3` | `-T cichlid_tree`: the same IQ-TREE command plus a fixed seed and crash recovery (**pending approval**) |
| `-s <tree>` | Use our species tree | the catalog's tree, unedited | same, checked by checksum |
| `-t <n>` | Parallel jobs | the machine's business | the core count (32) |
| `-A famsa` | Aligner | OrthoFinder's default | default (FAMSA), not changed |
| everything else | | OrthoFinder's defaults | defaults |

**IQ-TREE**, once per family (built by `cichlid-orthofinder-tree`; the template is recorded in `run_settings.json`):

```
iqtree3 -s <alignment> --prefix <name> -m MFP -seed 1 -T <threads> -quiet
```

| Flag | Meaning | The sheet | Ours |
|---|---|---|---|
| `-m MFP` | ModelFinder chooses the best model for each family | ModelFinder on every cluster | `-m MFP` for every family (default) |
| `-seed 1` | Fixed random seed, so a rerun gives the same tree | a fixed seed, recorded (already agreed) | `-seed 1`, recorded |
| `-T <n>` | Threads for this tree | the machine's business | 1, or more for large families (**pending approval**, see below) |
| `--prefix`, `-quiet` | Output file names, less screen output | — | no effect on the tree |

### IQ-TREE model options

| `-m` setting | What it does | Effect |
|---|---|---|
| **`MFP`** (the sheet's setting; our default) | Tests the full set of models (1,232 in our check) and keeps the best | best fit; the slowest |
| `TEST` | Tests only the common models (224 in our check); skips the slow, flexible "FreeRate" models | faster for large families; close to best, slightly lower for large ones |
| A fixed model, e.g. `LG+G4` | No model choice | faster; lower accuracy (in the pilot, it was the best fit for 1 family in 25) |
| `--fast` (tree search) | Stops searching for a better tree sooner | faster; some loss on large families |
| FastTree (`-T fasttree`) | A different, approximate program | much faster; clearly lower on large families; two methods in one data set |

**Our recommendation (not final; Coordinator approval needed):** keep `-m MFP` for every family except the 19
largest (over 1,000 genes each, 0.1% of families), which would get `-m TEST` (`--model hybrid`, rule `MFP;1000:TEST`).
These 19 set the end of the run. The accuracy cost is likely small and limited to those 19; it has not been measured.
The choice can be changed while the run goes on, before those families start (`cichlid-orthofinder set-model`).
**Until it is approved, the run uses `-m MFP` for all.**

### Pending approval: the run-time features and their options

None of these changes the inputs or the output format. Each one changes only how the run is carried out.

| Feature | Option and default | What it does | Effect on the results | Approval |
|---|---|---|---|---|
| Tree program | `-T cichlid_tree` (OrthoFinder) | The same IQ-TREE command, plus the fixed seed, a store of finished trees, and resuming from IQ-TREE's checkpoint files | none (tested: same trees; a crash loses no finished tree) | Coordinator approval needed |
| Fixed seed | `--seed 1` (default) | Same seed for every tree, recorded | makes reruns reproducible | already agreed by the coordinator |
| Model by family size | `--model full` (default) / `hybrid` / a rule such as `MFP;1000:TEST`; changeable during the run with `set-model` | Which `-m` each family gets | none with `full`; slightly lower accuracy for the 19 largest with `hybrid` | `full`: the sheet's setting. `hybrid`: Coordinator approval needed |
| Threads per tree | `--thread-rule "1;600:2;1500:8;3000:16"` (default); `--cores` caps the total (default: all cores) | 1 thread per tree, 2 from 600 genes, 8 from 1,500, 16 from 3,000; never more than the cap | **Trees built with more than 1 thread can differ slightly between runs, even with the fixed seed** (tested). With `--thread-rule 1`, every tree is exactly reproducible but the largest are much slower | Coordinator approval needed |
| Job order | `--order ascending` (default) / `descending` | Smallest families first, using OrthoFinder's own option (a one-word fix in memory) | none (tested: byte-identical output) | Coordinator approval needed |
| Early results | `--snapshot-sizes 300,1000` (default) / `''` for none | The results, in the coordinator's format, for every family up to 300 and 1,000 genes, while the run continues (kept locally) | none (made from a copy; matched the final results in the test) | Coordinator approval needed (where they are shared) |
| Accepting a run | always on | A run is accepted only if every family of 3+ genes got its tree | none; it only catches failures | — |
| Family settings | `convert --tree-set raw_pruned` (default) / `resolved`; `--min-size 4` (default) | Which tree goes in each family folder; the smallest family that gets a folder | Orthogroups of 2–3 genes get no folder (they are listed) | Coordinator approval needed |
| IQ-TREE version | 3.1.3, as the sheet names (3.1.4 is now out) | — | moving would rebuild every finished tree | stay, unless the coordinator asks |

---

## Where to find what

| You want | Look in |
|---|---|
| The coordinator's instructions for this repository | `docs/orthofinder_repo_spec.md` (a verbatim copy; the original is in the coordinator's `docs/`) |
| The commands for the run, from checks to hand-over | `docs/RUNBOOK.md` |
| What the pilots found (OrthoFinder's folders, labels, timings, tests) | `docs/pilot.md` |
| Why and how OrthoFinder calls our tree program | `docs/TREE_PROGRAM.md` |
| The sheet's acceptance criteria, one by one | `docs/COMPLIANCE.md` |
| Which coordinator and sheet versions this follows | `docs/SPEC_VERSION.md` |
| What changed, and why | `CHANGES.md` |
| Starting, resuming and recording the run | `src/cichlid_orthofinder/run.py` |
| Converting the results to the coordinator's format | `src/cichlid_orthofinder/convert.py` |
| The copy of OrthoFinder's folder for `native/` | `src/cichlid_orthofinder/native.py` |
| The checks before a restart | `src/cichlid_orthofinder/verify.py` |
| How each gene tree is built and kept | `src/cichlid_orthofinder/treecache.py` |
| Smallest-first order | `src/cichlid_orthofinder/of_launcher.py` |
| Early results | `src/cichlid_orthofinder/snapshot.py`; converted by `convert.py` (`convert_snapshot`) |
| The progress report | `src/cichlid_orthofinder/status.py` |
| The commands (`cichlid-orthofinder run`, `convert`, `native`, `verify`, `status`, `set-model`, ...) | `src/cichlid_orthofinder/cli.py` |
| The tests, and a small real OrthoFinder output used by them | `tests/`, `tests/fixtures/pilot2/` |

---

## Install, run, test

```bash
conda activate phylo-trees                      # OrthoFinder 3.1.5, IQ-TREE 3.1.3 and the other programs
python -m pip install -e <coordinator clone>    # cichlid-coordinator>=0.14,<0.15
python -m pip install -e .

cichlid-orthofinder verify --work <work folder>                     # every gate must PASS
cichlid-orthofinder run --work <work folder> --resume --dry-run     # shows the command; starts nothing
cichlid-orthofinder run --work <work folder> --resume --detach      # the restart, once approved
cichlid-orthofinder status <work folder>/v1/our_species_tree        # progress

PYTHONPATH=src python -m unittest discover tests                    # the offline tests
```

The full sequence, including conversion, validation and hand-over, is in `docs/RUNBOOK.md`.
