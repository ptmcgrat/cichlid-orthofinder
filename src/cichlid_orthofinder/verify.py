"""The gates to pass before the full run restarts, each printed as PASS or FAIL and saved to gates.txt.

    cichlid-orthofinder verify --work <scratch> [--model full|hybrid|<rule>] [--out DIR]

1. the model setting in effect is valid: the rule is printed (and whether it departs from the sheet's,
   ModelFinder on every cluster); PASS if it parses;
2. `cichlid-data check --release v1 --part OrthoFinder --fetch` passes;
3. the catalog's proteomes_md5 of the part's species is the one the saved search was made from;
4. the data root is the real one: $CICHLID_DATA_ROOT as set (or the coordinator's configured root), and
   `cichlid-data where proteome Mzeb` inside it, not a pilot or scratch root;
5. every program version is readable, and the free disk is shown (the run needs roughly 20 to 30 GB more);
6. `cichlid-orthofinder run ... --resume --dry-run` runs; its output is saved, to send to the reviewer;
7. the repository is clean, and its commit is tagged and pushed, so the commit recorded is the code that ran;
8. the model settings: IQ-TREE on one small alignment with -m MFP and with -m TEST, reporting the
   "ModelFinder will test up to N protein models" line and the FreeRate (+R) models each one fitted;
9. every local proteome, re-hashed here (hashlib), has the MD5 the catalog records for that file;
10. the catalog's species for the part, the part's own species list in the catalog, and the part's list in the
    coordinator's config/releases.toml are the same 35, in the same order (the difference is printed);
11. the species tree the run gives -s has the MD5 the catalog recorded for that file when it was built (and
    the run's recorded tree text is the same file).

Gate 6 needs the reviewer: this file records and checks, it does not decide.
Exit code 0 only if every gate passes.
"""

import argparse
import gzip
import hashlib
import json
import io
import random
import shutil
import subprocess
import sys
import tempfile
import tomllib
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path

RELEASE, PART, RUN_NAME = "v1", "OrthoFinder", "our_species_tree"
EXPECTED_PROTEOMES_MD5 = "5df42f6efb77659a03c6a82d6d9c7b26"
NEEDED_GB = 30
SCRATCH_MARKERS = ("pilot", "scratch", "/tmp/")
# keys of IQ-TREE's model checkpoint (<prefix>.model.gz) that are not models
CHECKPOINT_KEYS = {"finishedFastMLTree", "initTree", "finishedModelInit", "finishedCandidateSet",
                   "finishedModelFinal", "best_model_AIC", "best_model_AICc", "best_model_BIC"}


class Gate:
    def __init__(self, number, name):
        self.number, self.name, self.ok, self.lines = number, name, False, []

    def note(self, line):
        self.lines.append(line)

    def show(self, say=print):
        say(f"gate {self.number} {'PASS' if self.ok else 'FAIL'}  {self.name}")
        for line in self.lines:
            say(f"         {line}")


def _run(cmd, timeout=1800):
    """-> (exit code, stdout + stderr) of a command, or (None, the error) if it could not start."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, str(e)


def gate_model(model):
    from .run import describe_model, model_rule
    g = Gate(1, "the model setting in effect is valid")
    try:
        g.note(describe_model(model_rule(model)))
        g.ok = True
    except ValueError as e:
        g.note(f"{model!r}: {e}")
    return g


def gate_check(out_dir):
    g = Gate(2, f"cichlid-data check --release {RELEASE} --part {PART} --fetch")
    rc, text = _run(["cichlid-data", "check", "--release", RELEASE, "--part", PART, "--fetch"])
    (out_dir / "cichlid_data_check.txt").write_text(text)
    g.ok = rc == 0
    g.note(f"exit {rc}; last line: {text.strip().splitlines()[-1] if text.strip() else '(no output)'}")
    return g


def gate_proteomes():
    from cichlid_data import Catalog, get_layout
    g = Gate(3, "the proteomes are the ones the saved search was made from")
    catalog = Catalog.load(get_layout())
    species = [s.prefix for s in catalog.species(release=RELEASE, part=PART)]
    md5 = catalog.proteomes_md5(species)
    g.ok = md5 == EXPECTED_PROTEOMES_MD5
    g.note(f"{len(species)} species; proteomes_md5 {md5} (expected {EXPECTED_PROTEOMES_MD5})")
    return g


def gate_data_root(environ):
    from cichlid_data import get_layout
    g = Gate(4, "the data root is the real one, not a pilot or scratch root")
    env_root = environ.get("CICHLID_DATA_ROOT", "")
    root = Path(get_layout().local_root).expanduser().resolve()
    rc, text = _run(["cichlid-data", "where", "proteome", "Mzeb"], timeout=120)
    where = Path(text.strip().splitlines()[-1]).resolve() if rc == 0 and text.strip() else None
    g.note(f"$CICHLID_DATA_ROOT = {env_root or '(not set: the coordinator uses its configured root)'}")
    g.note(f"coordinator's root: {root}")
    g.note(f"cichlid-data where proteome Mzeb: {where if where else text.strip() or '(no output)'}")
    inside = where is not None and root in where.parents
    scratch = [m for m in SCRATCH_MARKERS if m in str(root) or (where and m in str(where))]
    env_ok = not env_root or Path(env_root).expanduser().resolve() == root
    if not inside:
        g.note("the proteome is not inside the coordinator's root")
    if scratch:
        g.note(f"looks like a pilot or scratch root (contains {', '.join(scratch)})")
    if not env_ok:
        g.note("$CICHLID_DATA_ROOT differs from the coordinator's root")
    g.ok = inside and not scratch and env_ok and where is not None and where.is_file()
    return g


def gate_versions(work):
    from .run import missing_versions, program_versions
    g = Gate(5, "every program version is readable; free disk")
    versions = program_versions()
    g.note(", ".join(f"{k} {v or 'MISSING'}" for k, v in versions.items()))
    free_gb = shutil.disk_usage(Path(work).expanduser().resolve()).free / 1e9
    g.note(f"free disk at the work folder: {free_gb:.0f} GB (the run needs roughly 20 to {NEEDED_GB} GB more)")
    g.ok = not missing_versions(versions) and free_gb >= NEEDED_GB
    return g


def gate_dry_run(work, model, out_dir):
    from .run import run
    g = Gate(6, "cichlid-orthofinder run --resume --dry-run (send the output to the reviewer)")
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            run(work, run_name=RUN_NAME, release=RELEASE, model=model, resume=True, dry_run=True)
        g.ok = True
    except Exception as e:                       # report any refusal as the gate's result
        buf.write(f"\nREFUSED: {type(e).__name__}: {e}\n")
    (out_dir / "dry_run.txt").write_text(buf.getvalue())
    command = [l for l in buf.getvalue().splitlines() if l.startswith("would run:")]
    g.note(command[0] if command else buf.getvalue().strip().splitlines()[-1] if buf.getvalue().strip() else "")
    g.note(f"saved to {out_dir / 'dry_run.txt'}")
    return g


def gate_repo():
    g = Gate(7, "the repository is clean, and its commit is tagged and pushed")
    repo = str(Path(__file__).resolve().parents[2])
    rc, _ = _run(["git", "-C", repo, "rev-parse", "--is-inside-work-tree"], timeout=60)
    if rc != 0:
        g.note(f"{repo} is not a git repository: no commit, tag or push to record")
        g.note("expected to FAIL for now: git comes later, by the owner's instruction")
        return g
    _, head = _run(["git", "-C", repo, "rev-parse", "HEAD"], timeout=60)
    _, status = _run(["git", "-C", repo, "status", "--porcelain"], timeout=60)
    _, tags = _run(["git", "-C", repo, "tag", "--points-at", "HEAD"], timeout=60)
    _, remote = _run(["git", "-C", repo, "branch", "-r", "--contains", "HEAD"], timeout=60)
    clean, tagged, pushed = not status.strip(), bool(tags.strip()), bool(remote.strip())
    g.note(f"commit {head.strip()}")
    g.note(f"clean: {'yes' if clean else 'no (' + str(len(status.splitlines())) + ' changed)'}; "
           f"tag: {tags.strip() or 'none'}; on a remote branch: {remote.strip() or 'no'}")
    g.ok = clean and tagged and pushed
    return g


def small_alignment(n=6, length=120, seed=1):
    """A small protein alignment, the same every time: one sequence and n-1 copies with changes."""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    rng = random.Random(seed)
    base = [rng.choice(aa) for _ in range(length)]
    rows = []
    for i in range(n):
        s = list(base)
        for _ in range(8 * (i + 1)):
            s[rng.randrange(length)] = rng.choice(aa)
        rows.append(f">s{i}\n{''.join(s)}\n")
    return "".join(rows)


def gate_models(iqtree="iqtree3"):
    g = Gate(8, "model settings: what -m MFP and -m TEST test")
    found = {}
    with tempfile.TemporaryDirectory() as d:
        aln = Path(d) / "aln.fa"
        aln.write_text(small_alignment())
        for model in ("MFP", "TEST"):
            prefix = Path(d) / model
            rc, _ = _run([iqtree, "-s", str(aln), "-m", model, "-seed", "1", "-T", "1", "--prefix", str(prefix),
                          "-quiet"], timeout=600)
            log = Path(f"{prefix}.log")
            text = log.read_text(errors="replace") if log.is_file() else ""
            line = next((l.strip() for l in text.splitlines() if "ModelFinder will test up to" in l), None)
            fitted, freerate = [], 0
            if Path(f"{prefix}.model.gz").is_file():
                with gzip.open(f"{prefix}.model.gz", "rt") as fh:
                    fitted = [k for k in (l.split(":", 1)[0].strip() for l in fh if ":" in l)
                              if k[:1].isalpha() and "_" not in k and k not in CHECKPOINT_KEYS]
                freerate = sum("+R" in m for m in fitted)
            best = next((l.split(":", 1)[1].strip() for l in text.splitlines() if l.startswith("Best-fit model")), "?")
            found[model] = line
            g.note(f"-m {model}: {line or 'no ModelFinder line (exit ' + str(rc) + ')'}; "
                   f"{len(fitted)} models fitted, {freerate} FreeRate (+R); best {best}")
    g.ok = all(found.values())
    return g


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gate_proteome_files():
    """Each proteome re-hashed here, against the catalog's per-file proteome_md5."""
    from cichlid_data import Catalog, get_layout
    g = Gate(9, "every proteome file, re-hashed, has the catalog's MD5 for it")
    layout = get_layout()
    catalog = Catalog.load(layout)
    species = [s.prefix for s in catalog.species(release=RELEASE, part=PART)]
    bad = []
    for p in species:
        path = layout.proteome(p).local
        got = _md5(path) if Path(path).is_file() else "missing"
        if got != catalog[p].proteome_md5:
            bad.append(f"{p}: file {got}, catalog {catalog[p].proteome_md5}")
    g.note(f"{len(species)} proteomes re-hashed; mismatches: {len(bad)}")
    for b in bad:
        g.note(b)
    g.ok = bool(species) and not bad
    return g


def gate_species_lists():
    """The part's species three ways: the catalog's organisms for it, the part's list in the catalog, and the
    part's list in the coordinator's config/releases.toml (beside the part's species tree)."""
    from cichlid_data import Catalog, get_layout
    g = Gate(10, "the 35 species: catalog, the part's list, and the coordinator's releases.toml agree")
    catalog = Catalog.load(get_layout())
    in_catalog = [s.prefix for s in catalog.species(release=RELEASE, part=PART)]
    part_list = list(catalog.part(RELEASE, PART)["species"])
    releases_toml = Path(catalog.species_tree(RELEASE, PART)).parent / "releases.toml"
    config_list = tomllib.loads(releases_toml.read_text())["releases"][RELEASE]["parts"][PART]["species"]
    g.note(f"catalog species for the part: {len(in_catalog)}; the part's list: {len(part_list)}; "
           f"{releases_toml}: {len(config_list)}")
    for name, other in (("the part's list in the catalog", part_list), ("config/releases.toml", config_list)):
        only_cat = [p for p in in_catalog if p not in other]
        only_other = [p for p in other if p not in in_catalog]
        order = "same order" if other == in_catalog else "different order" if not (only_cat or only_other) else ""
        g.note(f"diff with {name}: only in the catalog {only_cat or 'none'}; only in {name} {only_other or 'none'}"
               + (f"; {order}" if order else ""))
    g.ok = len(in_catalog) == 35 and in_catalog == part_list == list(config_list)
    return g


def gate_species_tree(work):
    """The -s tree the run uses (run_settings.json, last attempt), against the checksum the catalog recorded
    for that file when it was built, and against the tree text the run recorded when it started."""
    from cichlid_data import Catalog, get_layout
    g = Gate(11, "the species tree used for -s is the catalog's, by checksum")
    catalog = Catalog.load(get_layout())
    name = catalog.part(RELEASE, PART).get("species_tree")
    recorded = catalog.meta.get("config", {}).get(name)
    settings_file = Path(work).expanduser().resolve() / RELEASE / RUN_NAME / "run_settings.json"
    if not settings_file.is_file():
        g.note(f"no run_settings.json at {settings_file}")
        return g
    s = json.loads(settings_file.read_text())
    used = Path(s["attempts"][-1]["settings"]["species_tree"])
    used_md5 = _md5(used) if used.is_file() else None
    text = s.get("species_tree_text")
    text_md5 = hashlib.md5((text + "\n").encode()).hexdigest() if text else None
    g.note(f"catalog's recorded MD5 for {name}: {recorded}")
    g.note(f"-s file of the run: {used} -> {used_md5 or 'missing'}")
    g.note(f"the run's recorded tree text (+ newline): {text_md5}")
    g.note(f"the catalog's tree now: {catalog.species_tree(RELEASE, PART)} -> {_md5(catalog.species_tree(RELEASE, PART))}")
    g.ok = recorded is not None and used_md5 == recorded and text_md5 == recorded
    return g


def main(argv=None):
    p = argparse.ArgumentParser(prog="cichlid-orthofinder verify", description=__doc__.split("\n\n")[0])
    p.add_argument("--work", required=True, help="the run's scratch folder (as given to cichlid-orthofinder run)")
    p.add_argument("--model", default="full",
                   help="the model setting to restart with: full (the sheet's setting; default), hybrid or a rule")
    p.add_argument("--out", help="folder for the saved outputs (default: a new temporary folder)")
    args = p.parse_args(argv)
    import os
    out_dir = Path(args.out or tempfile.mkdtemp(prefix="verify_")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    gates = [gate_model(args.model), gate_check(out_dir), gate_proteomes(),
             gate_data_root(os.environ), gate_versions(args.work), gate_dry_run(args.work, args.model, out_dir),
             gate_repo(), gate_models(), gate_proteome_files(), gate_species_lists(), gate_species_tree(args.work)]
    report = io.StringIO()
    for g in gates:
        g.show(lambda s: (print(s), report.write(s + "\n")))
    passed = sum(g.ok for g in gates)
    summary = f"{passed} of {len(gates)} gates pass" + ("" if passed == len(gates) else ": do not restart")
    print(summary)
    report.write(summary + "\n")
    (out_dir / "gates.txt").write_text(report.getvalue())
    print(f"outputs saved in {out_dir}")
    return 0 if passed == len(gates) else 1


if __name__ == "__main__":
    sys.exit(main())
