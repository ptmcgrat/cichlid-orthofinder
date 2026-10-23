# OrthoFinder repository: instruction sheet

*For the `cichlid-orthofinder` repository. October 2026. Written against cichlid-coordinator v0.13.0. Read first: `project_manifest.md` (why the project exists), `species_final.md` (the species), and `coordinator_spec.md` (the layout, the schemas and the family contract, which this sheet relies on).*

## 1. Role

This repository runs OrthoFinder on the project's 35 proteomes and turns its gene trees into the project's standard family format. It owns:

1. **The OrthoFinder run:** the settings in section 4, run reproducibly.
2. **The conversion:** OrthoFinder's orthogroups and gene trees into `Results/v1/OrthoFinder/assign/` and `families/`, following the family contract.
3. **OrthoFinder's own outputs,** kept as written in `Results/v1/OrthoFinder/native/`, so they are backed up and can be pulled by anyone.

It does **not** make orthology calls, name genes, find human or zebrafish orthologs, or analyze results. The homology decider (project 3, `cichlid-homology`) reads the families this repository writes and calls them, with the same rule it uses for TreeFam and Compara, so the three methods' calls are comparable: differences between them then reflect the trees, not the calling logic.

**Rules (the same as for every repository):**

- Never build a data path by hand. Ask the coordinator (`get_layout()`, `cichlid-data where ...`).
- Never call rclone. Moving data to and from the cloud is the coordinator's (`get_store()`, `cichlid-data pull | push | ensure`).
- Never download an input. The proteomes come from the catalog.
- Never copy the species list or the species tree into this repository. Read them from the catalog.

Depend on the coordinator with a range pinned to one minor version, for example `cichlid-coordinator>=0.13,<0.14` (the minor you have tested), and use only its public API (`cichlid_data`) and its command line. Below 1.0.0 a minor version may change the layout API (v0.9.0 did: results are now `Results/<release>/<part>/`). When the coordinator moves to a new minor version, run your tests and `cichlid-data validate` against it, then move your pin to it; `cichlid_data.__version__` is the installed version.

## 2. Repository layout (suggested)

```
cichlid-orthofinder/
├── pyproject.toml           package cichlid_orthofinder; depends on cichlid-coordinator>=0.13,<0.14
├── src/cichlid_orthofinder/
│   ├── run.py               prepare the input folder, run OrthoFinder, record the settings
│   ├── convert.py           orthogroups and gene trees -> assign/ and families/
│   ├── native.py            what is kept under native/
│   └── cli.py               the cichlid-orthofinder command
├── tests/                   offline: small fake OrthoFinder output folders
└── docs/                    this sheet, notes on the runs
```

Moving code in from elsewhere follows the same rule as in the coordinator: verbatim, with its tests, changing only paths, transfers and formats.

## 3. Inputs (all through the coordinator)

| Input | How to get it |
|---|---|
| The species | `Catalog.load(get_layout()).species(release="v1", part="OrthoFinder")`: the 35 species, in the catalog's order |
| The proteomes | `get_layout().proteome(prefix)`: one protein per gene, headers `>gene_id`, already valid. Make sure they are local first: `cichlid-data check --release v1 --part OrthoFinder --fetch` (which also verifies the species tree) |
| The species tree | `Catalog.load(layout).species_tree("v1", "OrthoFinder")`, or `cichlid-data catalog tree --release v1 --part OrthoFinder`: rooted, fully resolved, leaves = the species prefixes |

Run `cichlid-data check --release v1 --part OrthoFinder --fetch` before every run; stop if it fails.

**Names matter.** OrthoFinder takes a species' name from its proteome file name, and `-s` needs the species tree's leaves to be those names. The coordinator's proteomes are already `<prefix>.faa`. Do not rename a species anywhere.

**Inputs are read-only: use a scratch input folder and a separate output folder.** By default OrthoFinder writes its results inside the folder it reads. Pointed at the coordinator's `Inputs/Proteins/proteomes/`, it would write into an input folder (and read the proteomes of every organism there, not just the 35). So build a scratch input folder holding links to exactly the 35 proteomes from the catalog (`layout.proteome(prefix)`), and send the results to a scratch output folder (OrthoFinder's `-o`, a parent folder: the results go to `<-o>/Results_<date>/`; verified in the pilot). When the run is finished, copy the result folder into `native/` (section 5). Never write into `Inputs/`.

**Gene IDs.** Our gene IDs (`Mzeb_101467881`) contain an underscore, which is fine for the coordinator. OrthoFinder renames every sequence internally (`<species number>_<sequence number>`) and writes the mapping back to our names in its working directory (`SequenceIDs.txt` and `SpeciesIDs.txt`; verified in the pilot). The conversion in section 6 maps every name back with it. **OrthoFinder does put the species name in front of the ID** in the resolved trees and the alignments, so that `Mzeb_101467881` becomes `Mzeb_Mzeb_101467881` (`-X` stops that, but it is not an agreed setting; see "What the installed OrthoFinder does"). Whatever it does, the converted leaf must be exactly `<prefix>_<id>`.

## 4. The OrthoFinder runs

Agreed settings:

| Setting | Value |
|---|---|
| Gene trees | multiple sequence alignments and a maximum-likelihood tree for every cluster: `-M msa` (the default in 3.1.5) with IQ-TREE: **`-T iqtree3`** (verified in the pilot) |
| Species tree | our species tree, with `-s <species tree>` (the file from the catalog; do not edit it) |
| Second run | **without** `-s`, so OrthoFinder infers its own species tree; kept for comparison only. The owner decides whether and how it is run |
| Random seed | **a fixed `-seed`**, recorded in `run.json` (`parameters`): without it IQ-TREE draws a new seed on every run and the trees differ (decided by the coordinator, for reproducibility) |
| Model | OrthoFinder's built-in IQ-TREE command gives no model, so **every cluster runs ModelFinder**. That stays the setting; any other model setting needs the owner's approval and is recorded in `run.json` |
| Everything else | OrthoFinder's defaults, unless agreed otherwise. Thread counts are the machine's business and are not science |

**A pilot comes first.** Before the 35-species run, do a pilot on a handful of species (or subsampled proteomes) to see what OrthoFinder really writes: the folder structure, the mapping files, the labels in the trees, the option names of the installed version. Write the conversion and its tests against that real output. A 35-species run with `-M msa` and IQ-TREE on every orthogroup is heavy; **before it starts, estimate its run time, memory and disk from the pilot, say where it will run, and get the owner's agreement.**

### What the installed OrthoFinder does (verified in the pilot, 2026-10-21)

**Versions:** OrthoFinder 3.1.5, IQ-TREE 3.1.3, FAMSA 2.4.1 (the default aligner), DIAMOND 2.2.8, MCL 22-282.

**Options.** `-T` takes `fasttree`, `fasttree_fastest`, `raxml` or `iqtree3`. `-o` is a **parent** folder: results go to `<-o>/Results_<Mon><day>/`. `-X` stops species names being added to sequence IDs (not used: it is not an agreed setting). Undocumented but useful: `-og` (stop after the orthogroups), `-ogt` (stop after the trees), `-fg` (start from orthogroups), `-fgt` (start from gene trees), `-mt` (threads per external command; `-mt 4` with `-t 16` ran 64 threads, so it does not limit the total). The built-in IQ-TREE command is `iqtree3 -s INPUT --prefix PATH/IDENTIFIER -quiet`: no model (so ModelFinder on every tree, single-threaded) and a new random seed each run (a rerun gave the same clusters and byte-identical alignments, but 108 of 129 trees differed).

**Folders** of `Results_<date>/`: `Orthogroups/`, `Orthogroup_Sequences/`, `MultipleSequenceAlignments/`, `Resolved_Gene_Trees/` (**one file**, `Resolved_Gene_Trees.txt`, one line per orthogroup), `Species_Tree/`, `Phylogenetic_Hierarchical_Orthogroups/` (`N1.tsv` to `N3.tsv`; N0 only as `WorkingDirectory/N0.ids.tsv`), `Gene_Duplication_Events/`, `Orthologues/`, `Putative_Xenologs/`, `Comparative_Genomics_Statistics/`, `WorkingDirectory/`. **There is no `Gene_Trees/` folder.** `WorkingDirectory/` holds the search inputs and results (`Species<i>.fa`, `diamondDBSpecies<i>.dmnd`, `Blast<i>_<j>.txt.gz`), `clusters_OrthoFinder_I1.2.txt`, `Sequences_ids/`, `Alignments_ids/` (the alignments and IQ-TREE's own files) and `Trees_ids/` (the trees as inferred), `pickle/`, `N0.ids.tsv`, and the mapping files.

**Mapping files** (`WorkingDirectory/`): `SpeciesIDs.txt` (`0: Acal.faa`, numbered **alphabetically**, not in catalog order) and `SequenceIDs.txt` (`0_0: Acal_113014445`, our IDs kept whole).

**Labels:** `Orthogroups/*.txt` and `*.tsv` carry our ID unchanged (`Acal_113012997`); `Resolved_Gene_Trees.txt`, `MultipleSequenceAlignments/` and the duplication tables **add the species name** (`Acal_Acal_113012997`); `Trees_ids/`, `Alignments_ids/` and `N0.ids.tsv` use internal names (`0_393`). The conversion ends with exactly `<prefix>_<id>`.

**Two sets of trees, and they are not what the sheet assumed.** OrthoFinder 3.1.5 builds its IQ-TREE trees on **preliminary clusters**, then splits each into the final orthogroups using the rooted trees, and renumbers them. In pilot 2, raw tree `OG0000000` (20 genes) became the final orthogroups `OG0000012`, `OG0000065` and `OG0000107` (9, 6 and 5 genes); none of the 319 raw trees has the gene set of the final orthogroup with the same number, but every final orthogroup lies inside one raw tree.

| | `WorkingDirectory/Trees_ids/` (as inferred) | `Resolved_Gene_Trees/Resolved_Gene_Trees.txt` |
|---|---|---|
| One per | preliminary cluster | final orthogroup (matches `Orthogroups.txt`, 336 of 336) |
| Rooted | **no** (three children at the root) | **yes** (two children; species tree) |
| Labels | internal (`0_393`) | species name added (`Acal_Acal_…`) |
| Branch lengths, supports | lengths yes, **no supports** (no bootstrap) | lengths yes, no supports (internal nodes named `n0`, ...) |
| Topology | IQ-TREE's | **differs from the raw tree in 69 of 336 (21%)** |

**Which tree is the family tree: the coordinator's recommendation, for the owner to decide.** The trees as inferred (`Trees_ids/`) are per preliminary cluster, not per orthogroup, and the resolved trees (`Resolved_Gene_Trees.txt`) were rearranged by OrthoFinder with the species tree in 21% of cases, which would build OrthoFinder's own duplication logic into one side of the comparison of methods. **Recommendation: for each final orthogroup, the raw IQ-TREE tree pruned to that orthogroup's genes** (`rooted: false`, `tree_mode: denovo`; the homology decider roots it by reconciliation with the master species tree, like the other two methods). The resolved set stays in `native/`. Pruning always works. **A check to report first:** how many final orthogroups come from a cluster that was split, and, for those, how many species the pieces share. A split that separates genes of the same species separates paralogs, which the families would then no longer show together; if that is common, the owner should know before deciding. Record the choice in `run.json` (`parameters`).

Every setting actually used goes into `run.json` (`parameters`), together with the OrthoFinder version (`method_version`), the IQ-TREE and alignment-program versions, and the command line. A run that can't say how it was made is not accepted.

**Where the two runs go.** The first run, with our species tree, is the one that becomes `Results/v1/OrthoFinder/`: its families are converted and called. The second run, with OrthoFinder's own tree, is a comparison: keep its OrthoFinder output under `Results/v1/OrthoFinder/native/own_species_tree/`, and do not convert it. The first run's output goes beside it, in `native/our_species_tree/`. *(This arrangement is an assumption to confirm: if calls are wanted from the second run it becomes a run of its own, with its own name and catalog entry.)* Each top-level folder under `native/` is stored as a separate tar, so adding the second run does not repack the first. Loose files directly in `native/` go together into one extra tar that is repacked whenever any of them changes, so put files in subfolders.

## 5. Outputs

Everything below is under `get_layout().results("v1", "OrthoFinder")` (a `ResultsLayout`; `v1` is the release, which is the coordinator's to name, and `OrthoFinder` the part):

| Path | Content |
|---|---|
| `assign/<prefix>.assign.tsv` | one per species, the existing assignment schema: `gene_id` and `protein_length` for every gene; `best_family` is the gene's orthogroup, **empty for a gene in none**; the score, coverage and fusion columns **empty** (they belong to placement methods) |
| `families/<family>/` | one folder per orthogroup, the family contract below |
| `native/` | OrthoFinder's own result folders, as it wrote them (see below) |
| `run.json` | the run record: `release` (`v1`; it must match the folder it is in), `method` (`OrthoFinder`), `method_version`, `parameters`, `started`, `finished`, `coordinator_version`, `catalog_md5`, `proteomes_md5`, `species` (the 35 prefixes). Write it with `schemas.write_run`, taking `coordinator_version` from `cichlid_data.__version__`, `catalog_md5` from `Catalog.load(layout).md5()` and `proteomes_md5` from `Catalog.load(layout).proteomes_md5(species)` (never compute any of them by hand). `proteomes_md5` is the one that ties the run to exactly the proteomes it used; `catalog_md5` identifies one `catalog.json` file and changes whenever the catalog is rebuilt |

### The family contract (summary; the full text is in `coordinator_spec.md`)

For every **orthogroup** (OrthoFinder's primary orthogroups, `OG0001234`; they have gene trees), a folder named by its family ID (letters, digits, `.`, `_`, `-`) holding:

- **`tree.nwk`**: OrthoFinder's gene tree with **every leaf renamed back to our gene ID** (`<prefix>_<id>`). All leaves are ours: every sequence comes from the 35 proteomes, so there are no reference-protein IDs. A leaf that still looks like `0_123` is refused by validation. Keep the branch lengths. Neither tree set has support values, so there are none to keep.
- **`combined.aln`**: required. The orthogroup's alignment (aligned FASTA), rows renamed to the same gene IDs, with a row for every leaf. **Check every alignment:** in pilot 2 two 3-gene orthogroups had a one-row alignment named `None`; the conversion refuses an alignment that lacks a row for any leaf. (The name is the placement methods'; they write other files beside it, which OrthoFinder families do not have.)
- **`done.json`**: a small JSON object, written **last** (a folder without it is unfinished, and a rerun redoes it), holding `"release": "v1"`, `"rooted"`, `"tree_mode"` and the inputs hash that lets a rerun skip finished families. `"rooted"`: the raw trees are **unrooted** (three children at the root; pruning does not root them), so write `"rooted": false`. Write `"rooted": true` with `"rooting": "orthofinder species tree"` **only if the root of that `tree.nwk` really has two children** (validation checks it). If a tree arrives unrooted, write `"rooted": false`, and the homology decider roots it by reconciliation with the master species tree (not by midpoint). Never write `"rooted": true` to mean "probably". `"tree_mode"` says how the tree was made; the trees here are built from scratch from the alignment, so use `"denovo"`, the word the TreeFam and Compara parts use for trees built from scratch (`place` or `denovo`).

The families are the **orthogroups (`OG…`)**; the hierarchical orthogroups (`N0.HOG…`) are not converted and stay in `native/` with the rest of OrthoFinder's output. An orthogroup for which OrthoFinder wrote no gene tree gets no family folder; its genes still appear in `assign/` with their orthogroup, and in OrthoFinder's own tables. Report how many there are in the conversion's summary, so the gap is visible.

### What is kept under `native/`

`native/our_species_tree/` (and `native/own_species_tree/` for the comparison run) is **a curated copy of OrthoFinder's result folder**: everything the run wrote, except the search files of `WorkingDirectory/`. **Keep and push:** every top-level result folder (`Orthogroups/`, `Orthogroup_Sequences/`, `MultipleSequenceAlignments/`, `Resolved_Gene_Trees/`, `Species_Tree/`, `Phylogenetic_Hierarchical_Orthogroups/`, `Gene_Duplication_Events/`, `Orthologues/`, `Putative_Xenologs/`, `Comparative_Genomics_Statistics/`, the logs), and from `WorkingDirectory/` only `SpeciesIDs.txt`, `SequenceIDs.txt`, `N0.ids.tsv`, `clusters_OrthoFinder_I1.2.txt`, `Trees_ids/` and `Alignments_ids/` (the raw trees and alignments, which are hours of computation and are what the families are built from). **Do not push:** the search inputs and results (`Species<i>.fa`, `diamondDBSpecies<i>.dmnd`, `Blast<i>_<j>.txt.gz`: 8.7 GB for the 35 species, regenerable in about 3 h) or `pickle/`. Keep the search results on the machine that ran them, read-only, with their `proteomes_md5`; they are reusable with `-b`. You make the copy, so nothing needs to change in the coordinator. After the full run, report the size of each top-level folder; if one is too large to push, say so and the coordinator decides. Never push `native/` before the keep-list has been applied and the sizes reported.

## 6. Conversion (steps)

1. Read the mapping files for the run; check every species and sequence is accounted for (35 species; the count of sequences equals the count of proteins in the proteomes).
2. For every species, write `assign/<prefix>.assign.tsv` from the orthogroup table through `schemas.ASSIGN.write` (it refuses invalid rows). Every gene of the proteome appears, in the proteome's order.
3. For every orthogroup with a tree, write the family folder: renamed tree, renamed alignment, `done.json`. Write `done.json` **last**: a folder without it is unfinished, and a rerun redoes it.
4. Write `run.json` with `schemas.write_run`.
5. Validate (section 7). Then, and only then, push.

Conversion is deterministic: the same OrthoFinder output gives the same files, byte for byte.

## 7. Validation and handing over

```bash
cichlid-data check --release v1 --part OrthoFinder --fetch           # before a run: 35 proteomes and the species tree
R=$(cichlid-data where results v1 OrthoFinder)                        # never type the data root
cichlid-data validate "$R"/assign/*.assign.tsv                         # one per species
cichlid-data validate "$R"/families                                    # every family folder, leaf names checked against the catalog
cichlid-data validate "$R"/run.json
cichlid-data push 'results/v1/OrthoFinder/*' --dry-run                # then without --dry-run
```

Hand over to the homology decider only when all of these pass. The tests of this repository use small fake OrthoFinder output folders and run offline; they must cover the renaming (including a leaf that cannot be mapped, and the species-name prefix OrthoFinder may add), a gene in no orthogroup, an orthogroup with no gene tree, rooted and unrooted trees, and the byte-for-byte repeatability of the conversion. Build the fake folders from the pilot's real output, so they have its real structure. There are no earlier OrthoFinder results to compare with: the verification is the offline tests, `cichlid-data validate` on every output, and a `run.json` that records the version of every program.

## 8. What this repository must not do

- **No orthology calls.** No human or zebrafish relationships, no "one-to-one", no names. OrthoFinder's own ortholog tables stay in `native/`; they are not project results.
- **No hand-built paths,** and no copy of the species list or the species tree.
- **No rclone,** no downloads, no editing of the proteomes or of the catalog, and nothing written into an input folder.
- **No change to the species tree** to make a run work. If the tree needs changing, that is a decision recorded in the coordinator (`species_final.md` and `config/species_tree_orthofinder.nwk`, together).
- **No change to a file format.** If the contract is not enough, ask the coordinator; a schema change is a version change.

## 9. Acceptance criteria

1. The tests run offline and pass.
2. `cichlid-data check --release v1 --part OrthoFinder` passes, and the species tree used for `-s` is the catalog's, byte for byte.
3. `run.json` records the OrthoFinder version, every setting, the versions of IQ-TREE and the alignment program, the tree set chosen, and the coordinator version and catalog checksum, and passes `cichlid-data validate`.
4. Every `assign` file and every family folder passes `cichlid-data validate`, with no leaf left in OrthoFinder's internal naming.
5. Every gene of every proteome appears in its species' assignment file.
6. The first run's results and its native outputs are pushed, and `cichlid-data status 'results/v1/OrthoFinder/*'` shows everything in sync.
7. No other path, species list or cloud call appears anywhere in the code, and nothing is written into an input folder.
8. The pilot's findings (folder structure, mapping files, labels, option names) and the tree set chosen for the families are written down in the repository's `docs/` before the full run.
