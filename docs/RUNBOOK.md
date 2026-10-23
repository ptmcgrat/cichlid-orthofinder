# Runbook: the 35-species OrthoFinder run, from restart to hand-over

Every data path comes from the coordinator (`cichlid-data where ...`, `cichlid-data catalog tree ...`); none is
typed by hand. The only paths you choose are the run's scratch folder and the decisions in angle brackets.

| Placeholder | What it is |
|---|---|
| `<WORK>` | the run's scratch folder (the stopped run: `~/Desktop/orthofinder_claude/work`); not a data folder |
| `<RUN>` | `<WORK>/v1/our_species_tree`, the run's folder inside it |
| `<MODEL>` | the approved model setting: `full` (ModelFinder on every cluster, the sheet's setting) or `hybrid` |
| `<APPROVAL>` | the owner's written approval: `"<who>, <YYYY-MM-DD>, <what was approved>"` |
| `<TREE_SET>`, `<ALIGNMENT>`, `<MIN_SIZE>` | the owner's choices for the families (defaults `raw_pruned`, `orthogroup`, `4`) |

```bash
conda activate phylo-trees
```

## 1. Before the restart

```bash
cichlid-data check --release v1 --part OrthoFinder --fetch        # the 35 proteomes and the species tree
cichlid-orthofinder verify --work <WORK> --model <MODEL> --approval <APPROVAL> --out <WORK>/verify
#   every gate PASS (gates.txt in <WORK>/verify); gate 7 needs the code committed, tagged and pushed first
cichlid-orthofinder run --work <WORK> --model <MODEL> --resume --dry-run
#   prints the command and the environment; nothing is started. Send it to the reviewer.
```

## 2. The restart

```bash
cichlid-orthofinder run --work <WORK> --model <MODEL> --approval <APPROVAL> --resume --detach
```

`--resume` reuses the run's search and every finished tree in its cache. To start from scratch instead, use a
new `<WORK>` and add `--from-search <the saved search's WorkingDirectory>`: it is used only if it is complete and
of exactly these proteomes.

## 3. Watching it

```bash
cichlid-orthofinder status <RUN>                 # stage, work done, the three stage dates, the running trees
cichlid-orthofinder set-model <RUN> full         # only if the owner changes the model setting during the run
```

If the machine stops: the restart command of section 2 again (with `--resume`). Early results appear under
`<RUN>/snapshots/`; they are for looking, never converted or pushed.

## 4. After the run

The run is finished when `status` says `done` and the runner accepted it (`results` in `<RUN>/run_settings.json`).
`<RESULTS>` is that folder (also the first line of `status`).

```bash
R=$(cichlid-data where results v1 OrthoFinder)                   # never type the data root

cichlid-orthofinder convert --results <RESULTS> --run-settings <RUN>/run_settings.json --cache <RUN>/tree_cache \
    --tree-set <TREE_SET> --alignment <ALIGNMENT> --min-size <MIN_SIZE>
#   writes $R/assign/, $R/families/, $R/native/summary/ and $R/run.json, through the coordinator's layout

cichlid-orthofinder native --results <RESULTS> --dest "$(cichlid-data where results v1 OrthoFinder native)/our_species_tree"
#   prints the size of every top-level folder: report them before anything is pushed

cichlid-data validate "$R"/assign/*.assign.tsv
cichlid-data validate "$R"/families
cichlid-data validate "$R"/run.json

cichlid-data push 'results/v1/OrthoFinder/*' --dry-run           # review, then ask the owner
#   only with the owner's go-ahead:
cichlid-data push 'results/v1/OrthoFinder/*'
cichlid-data status 'results/v1/OrthoFinder/*'                   # everything in sync
```

The second run without `-s` (`--run own_species_tree`) only if the owner says yes; its output goes to
`native/own_species_tree/` with `cichlid-orthofinder native` and is not converted.
