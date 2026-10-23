# cichlid-orthofinder

Runs OrthoFinder on the project's 35 proteomes and converts its orthogroups and gene trees into the
project's family format (`Results/<release>/OrthoFinder/`). The instruction sheet is
`orthofinder_repo_spec.md` in the coordinator's `docs/`; pilot findings are in `docs/pilot.md`.

## Install

```bash
conda activate phylo-trees
python -m pip install -e <coordinator clone>    # cichlid-coordinator>=0.14,<0.15
python -m pip install -e .
```

## Tests

```bash
PYTHONPATH=src python -m unittest discover tests
```
