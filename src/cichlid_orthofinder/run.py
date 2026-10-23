"""Prepare the input folder, run OrthoFinder through the tree wrapper, record the settings.

Two runs, the same apart from the species tree (orthofinder_repo_spec.md, section 4):

    our_species_tree   -s <the catalog's species tree>; its families are converted
    own_species_tree   no -s, OrthoFinder infers its own species tree; comparison only

Everything comes from the coordinator: the species and their order from the catalog, the proteomes
from the layout, the species tree from the catalog. OrthoFinder reads a scratch folder of links to
exactly those proteomes and writes to a separate output folder, so nothing is written into Inputs/.

Gene trees go through the wrapper (treecache.py), registered in ~/config_orthofinder_user.json as the
tree program `cichlid_tree`: IQ-TREE with the agreed settings (ModelFinder), a fixed seed, threads by
cluster size within a core budget, and a cache in <work>/<release>/<run>/tree_cache. After a crash,
`--resume` restarts OrthoFinder from the finished search (-b); the cache returns every finished tree and
IQ-TREE resumes the unfinished ones from their checkpoints.

Layout of the work folder (scratch, given on the command line):
    <work>/<release>/<run>/input/          links to the proteomes
    <work>/<release>/<run>/output/         OrthoFinder's Results_<date>/ (a resume nests a new one inside)
    <work>/<release>/<run>/tree_cache/     the wrapper's cache and events.jsonl
    <work>/<release>/<run>/attempt_<n>.log OrthoFinder's output, per attempt
    <work>/<release>/<run>/run_settings.json
"""

import gzip
import hashlib
import json
import shutil
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .of_launcher import TREES_MSA_MD5

RUNS = ("our_species_tree", "own_species_tree")

# Agreed settings: -M msa, IQ-TREE (through the wrapper, which adds the fixed seed and the cache).
TREE_PROGRAM = "cichlid_tree"
METHOD = ["-M", "msa", "-T", TREE_PROGRAM]
TREE_ENTRY = {"program_type": "tree",
              "cmd_line": "cichlid-orthofinder-tree INPUT PATH/IDENTIFIER",
              "output_filename": "PATH/IDENTIFIER.treefile"}
USER_CONFIG = Path.home() / "config_orthofinder_user.json"

# program -> (command, pattern for the version in its output)
VERSIONS = {
    "orthofinder": (["orthofinder", "-v"], r"OrthoFinder:v?(\S+)"),
    "iqtree3": (["iqtree3", "--version"], r"IQ-TREE version (\S+)"),
    "famsa": (["famsa"], r"version (\S+)"),
    "diamond": (["diamond", "version"], r"diamond version (\S+)"),
    "mcl": (["mcl", "--version"], r"mcl (\S+)"),
}

ANSI = re.compile(r"\x1b\[[0-9;]*m")


def program_versions(runner=subprocess.run):
    """-> {program: version or None}, from each program's own report."""
    out = {}
    for name, (cmd, pattern) in VERSIONS.items():
        try:
            p = runner(cmd, capture_output=True, text=True, timeout=60)
            text = ANSI.sub("", (p.stdout or "") + (p.stderr or ""))
        except (OSError, subprocess.TimeoutExpired):
            text = ""
        m = re.search(pattern, text)
        out[name] = m.group(1) if m else None
    return out


def missing_versions(versions):
    """-> the programs whose version could not be read (the run refuses to start without all of them)."""
    return sorted(name for name, v in versions.items() if not v)


def repo_state(runner=subprocess.run):
    """-> {"commit": this repository's git commit or None, "dirty": uncommitted changes?}."""
    repo = Path(__file__).resolve().parents[2]
    try:
        head = runner(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=30)
        status = runner(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return {"commit": None, "dirty": None}
    commit = head.stdout.strip() if head.returncode == 0 else None
    return {"commit": commit or None, "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None}


def code_record(config=USER_CONFIG):
    """What ran: the package version, the SHA-256 of each of its source files, the tree program as registered
    with OrthoFinder and the IQ-TREE command the tree wrapper builds."""
    from . import __version__
    from .treecache import iqtree_args
    src = Path(__file__).resolve().parent
    hashes = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(src.glob("*.py"))}
    try:
        registered = json.loads(Path(config).read_text()).get(TREE_PROGRAM) if Path(config).is_file() else None
    except ValueError:
        registered = "unreadable"
    template = ["iqtree3", *iqtree_args("<alignment>", "<prefix>", "<seed>", "<threads>", "<model>")]
    return {"package_version": __version__, "source_sha256": hashes, "tree_program": TREE_ENTRY,
            "orthofinder_user_config_entry": registered, "iqtree_argument_template": template,
            "launcher_patch_md5": TREES_MSA_MD5, **repo_state()}


def register_tree_program(config=USER_CONFIG):
    """Make sure OrthoFinder's user configuration has the wrapper entry (other entries are kept)."""
    config = Path(config)
    data = json.loads(config.read_text()) if config.is_file() else {}
    if data.get(TREE_PROGRAM) != TREE_ENTRY:
        data[TREE_PROGRAM] = TREE_ENTRY
        config.write_text(json.dumps(data, indent=2) + "\n")
    return config


def prepare_input(folder, proteomes):
    """Link each species' proteome into `folder` as <prefix>.faa (OrthoFinder names a species by its
    file name). `proteomes` is an ordered {prefix: path}. The folder must be new or empty."""
    folder = Path(folder)
    if folder.exists() and any(folder.iterdir()):
        raise RuntimeError(f"input folder is not empty: {folder}")
    folder.mkdir(parents=True, exist_ok=True)
    for prefix, path in proteomes.items():
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"proteome for {prefix} not found: {path}")
        (folder / f"{prefix}.faa").symlink_to(path.resolve())
    return folder


LAUNCHER = "cichlid-orthofinder-of"     # OrthoFinder, with the tree order from $CICHLID_OF_ORDER (of_launcher.py)


def build_command(input_dir, output_dir, species_tree=None, threads=None, analysis_threads=None):
    """The OrthoFinder command line for a new run. `species_tree` None means OrthoFinder infers its own."""
    cmd = [LAUNCHER, "-f", str(input_dir), "-o", str(output_dir), *METHOD]
    return cmd + _common(species_tree, threads, analysis_threads)


def resume_command(working_dir, species_tree=None, threads=None, analysis_threads=None):
    """The command that restarts a run from its finished search (`-b`); the trees come from the cache."""
    return [LAUNCHER, "-b", str(working_dir), *METHOD] + _common(species_tree, threads, analysis_threads)


def _common(species_tree, threads, analysis_threads):
    cmd = []
    if species_tree is not None:
        cmd += ["-s", str(species_tree)]
    if threads:
        cmd += ["-t", str(threads)]
    if analysis_threads:
        cmd += ["-a", str(analysis_threads)]
    return cmd


def search_problems(working_dir):
    """What is missing or unfinished in a run's all-against-all search. -> [] when it is complete.

    The search writes one WorkingDirectory/Blast<i>_<j>.txt.gz per ordered pair of species (from
    SpeciesIDs.txt). A run killed during the search leaves some missing or cut short; resuming from such a
    search (`-b`) would cluster from partial hits."""
    wd = Path(working_dir)
    ids = wd / "SpeciesIDs.txt"
    if not ids.is_file():
        return ["SpeciesIDs.txt is missing (the search never started)"]
    n = sum(1 for l in ids.read_text().splitlines() if l.strip())
    problems = []
    for i in range(n):
        for j in range(n):
            f = wd / f"Blast{i}_{j}.txt.gz"
            if not f.is_file():
                problems.append(f"{f.name} missing")
                continue
            try:
                with gzip.open(f, "rb") as fh:
                    while fh.read(1 << 20):
                        pass
            except (OSError, EOFError) as e:
                problems.append(f"{f.name} unreadable ({e})")
    return problems


def _fasta(path):
    """-> [(first word of the header, sequence)] of a FASTA file, sequences upper-case, no final '*'."""
    out, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq).upper().rstrip("*")))
                name, seq = line[1:].split()[0] if line[1:].split() else "", []
            elif line:
                seq.append(line)
    if name is not None:
        out.append((name, "".join(seq).upper().rstrip("*")))
    return out


def search_mismatches(working_dir, proteomes):
    """Does a saved search belong to exactly these proteomes ({prefix: path})? -> [] when it does.

    Checks the species (SpeciesIDs.txt) and, species by species, that the search's renamed copy
    (Species<i>.fa) holds the same sequences, in the same order, as the proteome, and that SequenceIDs.txt
    maps each to the proteome's own ID."""
    wd = Path(working_dir)
    species = {}
    for l in (wd / "SpeciesIDs.txt").read_text().splitlines():
        if l.strip():
            i, name = l.split(":", 1)
            species[int(i)] = name.strip().removesuffix(".faa")
    if sorted(species.values()) != sorted(proteomes):
        return [f"species differ: search {sorted(species.values())}, proteomes {sorted(proteomes)}"]
    ids = {}
    for l in (wd / "SequenceIDs.txt").read_text().splitlines():
        if l.strip():
            k, v = l.split(":", 1)
            ids[k.strip()] = v.strip().split()[0]
    problems = []
    for i, prefix in sorted(species.items()):
        mine = _fasta(proteomes[prefix])
        theirs = _fasta(wd / f"Species{i}.fa")
        if len(mine) != len(theirs):
            problems.append(f"{prefix}: {len(theirs)} sequences in the search, {len(mine)} in the proteome")
            continue
        for j, ((pid, pseq), (sid, sseq)) in enumerate(zip(mine, theirs)):
            if sid != f"{i}_{j}" or ids.get(sid) != pid or sseq != pseq:
                problems.append(f"{prefix}: sequence {j + 1} ({pid}) differs from the search's {sid}")
                break
    return problems


SEARCH_FILES = ("SpeciesIDs.txt", "SequenceIDs.txt", "Species*.fa", "Blast*.txt.gz", "diamondDBSpecies*.dmnd")


def copy_search(src, dest):
    """Copy a finished search's files (not its clusters or trees) into a new WorkingDirectory. -> dest."""
    src, dest = Path(src), Path(dest)
    dest.mkdir(parents=True)
    for pattern in SEARCH_FILES:
        for f in sorted(src.glob(pattern)):
            shutil.copy2(f, dest / f.name)
    return dest


def set_aside(output_dir):
    """Move an unusable output folder out of the way (kept, not deleted). -> the new path."""
    aside = Path(output_dir).with_name(f"{Path(output_dir).name}_incomplete_{time.strftime('%Y%m%d_%H%M%S')}")
    Path(output_dir).rename(aside)
    return aside


def find_results(output_dir):
    """OrthoFinder writes one Results_<date> folder inside -o. -> that folder."""
    found = sorted(p for p in Path(output_dir).iterdir() if p.is_dir() and p.name.startswith("Results_"))
    if len(found) != 1:
        raise RuntimeError(f"expected one Results_* folder in {output_dir}, found {[p.name for p in found]}")
    return found[0]


def latest_results(output_dir):
    """The newest Results_* folder of a run, following resumes (each nests one in the previous
    WorkingDirectory/OrthoFinder/)."""
    r = find_results(output_dir)
    while True:
        nested = Path(r) / "WorkingDirectory" / "OrthoFinder"
        inner = sorted(nested.glob("Results_*")) if nested.is_dir() else []
        if not inner:
            return r
        r = inner[-1]


MODEL_SHORTCUTS = {"hybrid": "MFP;1000:TEST", "full": "MFP"}
# The sheet's agreed setting (section 4): ModelFinder on every cluster. Any other rule must be asked for.
DEFAULT_MODEL = "full"

# Attempts made before each attempt recorded its coordinator version (2026-10-23) ran under 0.13.0, the version
# installed on this machine until then; the run's first attempt (2026-10-22) is one of them.
COORDINATOR_RECORDED_FROM = "2026-10-23"
EARLIER_COORDINATOR = "0.13.0"


def note_earlier_coordinator(attempts):
    """Give every earlier attempt that did not record its coordinator version the one it ran under, with a note
    that it was recorded afterwards. Changes `attempts` in place; returns how many were filled in."""
    filled = 0
    for a in attempts:
        if "coordinator_version" not in a:
            earlier = str(a.get("started", "")) < COORDINATOR_RECORDED_FROM
            a["coordinator_version"] = EARLIER_COORDINATOR if earlier else None
            a["coordinator_version_note"] = (
                f"recorded afterwards: this attempt ran before attempts recorded their coordinator version "
                f"({COORDINATOR_RECORDED_FROM}), under cichlid-coordinator {EARLIER_COORDINATOR}" if earlier
                else "not recorded by this attempt")
            filled += 1
    return filled
LIVE_SETTINGS = "tree_settings.json"


def model_rule(rule):
    """'hybrid', 'full' or a rule such as 'MFP;1000:TEST' -> the checked rule."""
    from .treecache import parse_models
    rule = MODEL_SHORTCUTS.get(str(rule).strip().lower(), str(rule).strip())
    parse_models(rule)                       # raises ValueError if it is not a valid rule
    return rule


def describe_model(rule):
    """-> one line naming the model rule in effect, and saying so if it departs from the sheet's setting."""
    names = {v: k for k, v in MODEL_SHORTCUTS.items()}
    line = f"model setting: {rule}" + (f" ({names[rule]})" if rule in names else "")
    if rule != MODEL_SHORTCUTS[DEFAULT_MODEL]:
        line += ("; this departs from the sheet's setting, ModelFinder (-m MFP) on every cluster "
                 "(orthofinder_repo_spec.md, section 4)")
    return line


def read_model(work):
    """-> the model rule in a run's live settings file, or None."""
    f = Path(work) / LIVE_SETTINGS
    return json.loads(f.read_text())["model"] if f.is_file() else None


def write_model(work, rule):
    """Set the model rule that trees starting from now will use (atomic; keeps a history). -> the rule."""
    from .treecache import _atomic_write
    rule = model_rule(rule)
    f = Path(work) / LIVE_SETTINGS
    old = json.loads(f.read_text()) if f.is_file() else {}
    history = old.get("history", []) + ([{"model": old["model"], "since": old["since"]}] if "model" in old else [])
    _atomic_write(f, json.dumps({"model": rule, "since": _now(), "history": history}, indent=1) + "\n")
    return rule


def tree_env(cache, seed, thread_rule, cores, order="ascending", model="MFP", live=None):
    env = {"CICHLID_TREE_CACHE": str(cache), "CICHLID_TREE_SEED": str(seed),
           "CICHLID_TREE_THREADS": str(thread_rule), "CICHLID_TREE_CORES": str(cores),
           "CICHLID_OF_ORDER": order, "CICHLID_TREE_MODEL": model}
    if live:
        env["CICHLID_TREE_SETTINGS"] = str(live)
    return env


def _now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def run(work, run_name="our_species_tree", release="v1", threads=None, analysis_threads=None,
        seed=1, thread_rule="1;600:2;1500:8;3000:16", cores=None, order="ascending", snapshot_sizes=(300, 1000),
        model=DEFAULT_MODEL, from_search=None, poll_seconds=300,
        resume=False, dry_run=False, say=print):
    """Check the inputs, prepare the input folder (or resume), run OrthoFinder, record the attempt.
    -> the path of run_settings.json (or the command, for a dry run)."""
    import cichlid_data
    from cichlid_data import Catalog, get_layout, get_store

    if run_name not in RUNS:
        raise ValueError(f"run must be one of {RUNS}")
    model = model or DEFAULT_MODEL
    work = Path(work).expanduser().resolve() / release / run_name
    input_dir, output_dir, cache = work / "input", work / "output", work / "tree_cache"
    cores = int(cores or os.cpu_count() or 1)
    threads = threads or cores

    # 1. The coordinator's check: every proteome and the species tree. Stop if it fails.
    check = ["cichlid-data", "check", "--release", release, "--part", "OrthoFinder", "--fetch"]
    say("$ " + " ".join(check))
    if subprocess.run(check).returncode != 0:
        raise RuntimeError("cichlid-data check failed; not running")

    # 2. Species, proteomes and species tree, from the catalog (the cloud's, if there is none here).
    layout = get_layout()
    get_store().ensure(layout.catalog())
    catalog = Catalog.load(layout)
    species = [s.prefix for s in catalog.species(release=release, part="OrthoFinder")]
    proteomes = {p: layout.proteome(p).local for p in species}
    tree = catalog.species_tree(release, "OrthoFinder") if run_name == "our_species_tree" else None
    proteomes_md5 = catalog.proteomes_md5(species)

    settings_file = work / "run_settings.json"
    if resume:
        if not settings_file.is_file():
            raise RuntimeError(f"nothing to resume: {settings_file} does not exist")
        settings = json.loads(settings_file.read_text())
        note_earlier_coordinator(settings.get("attempts", []))   # written with the next attempt (not by a dry run)
        if settings["proteomes_md5"] != proteomes_md5:
            raise RuntimeError("the proteomes changed since this run started; start a new run instead")
        try:
            first = find_results(output_dir) / "WorkingDirectory"
            problems = search_problems(first)
        except (FileNotFoundError, RuntimeError) as e:   # no output folder, or no single Results_* in it
            first, problems = None, [str(e)]
        if problems:
            # killed during the search: a resume would cluster from partial hits, so start the search again
            say(f"the search did not finish ({len(problems)} problem(s), e.g. {problems[0]}); "
                "starting from the beginning; finished trees in the cache are still reused")
            restart_reason = {"search_problems": len(problems), "example": problems[0]}
            if output_dir.exists() and not dry_run:
                restart_reason["set_aside"] = str(set_aside(output_dir))
            cmd = build_command(input_dir, output_dir, tree, threads, analysis_threads)
        else:
            restart_reason = None
            cmd = resume_command(first, tree, threads, analysis_threads)
    else:
        if output_dir.exists():
            raise RuntimeError(f"output folder exists (a previous run? use --resume): {output_dir}")
        restart_reason = None
        reuse = None
        if from_search:
            # a finished search of exactly these proteomes: copied in, and OrthoFinder starts from it (-b)
            src = Path(from_search).expanduser().resolve()
            say(f"checking the saved search {src} (every file complete, the same proteomes)")
            problems = search_problems(src)
            if not problems:
                problems = search_mismatches(src, proteomes)
            if problems:
                # not complete, or not these proteomes: search again rather than reuse it
                say(f"the saved search is not reused ({len(problems)} problem(s), e.g. {problems[0]}); "
                    "searching again")
                reuse = {"from": str(src), "refused": problems[:20]}
            else:
                reuse = {"from": str(src), "to": str(output_dir / "Results_from_search" / "WorkingDirectory")}
        if reuse and "to" in reuse:
            cmd = resume_command(reuse["to"], tree, threads, analysis_threads)
        else:
            cmd = build_command(input_dir, output_dir, tree, threads, analysis_threads)
        settings = {
            "release": release, "run": run_name, "method": "OrthoFinder",
            "species": species, "catalog_md5": catalog.md5(), "proteomes_md5": proteomes_md5,
            "species_tree_text": Path(tree).read_text().strip() if tree else None,
            "search_reused": reuse,
            "attempts": [],
        }
    # The model rule: given, or the default (full). It is written to the live settings file, which trees read
    # again each time one starts, so `cichlid-orthofinder set-model` changes it while the run goes on.
    current = read_model(work) if resume else None
    rule = model_rule(model)
    say(describe_model(rule))
    env_add = tree_env(cache, seed, thread_rule, cores, order, rule, work / LIVE_SETTINGS)
    versions = program_versions()
    say("versions: " + ", ".join(f"{k} {v}" for k, v in versions.items()) + f"; coordinator {cichlid_data.__version__}")
    if missing_versions(versions):
        raise RuntimeError(f"could not read the version of {', '.join(missing_versions(versions))}; not running")
    code = code_record()
    if dry_run:
        say(f"{len(species)} species: {', '.join(species)}")
        say("environment: " + " ".join(f"{k}={v}" for k, v in env_add.items()))
        say("would run: " + " ".join(cmd))
        return cmd

    if not resume:
        prepare_input(input_dir, proteomes)
        if reuse and "to" in reuse:
            say(f"copying the saved search into {reuse['to']}")
            copy_search(reuse["from"], reuse["to"])
            problems = search_problems(reuse["to"])
            if problems:
                raise RuntimeError(f"the copied search is incomplete: {problems[0]}")
    elif restart_reason and not input_dir.is_dir():
        prepare_input(input_dir, proteomes)
    cache.mkdir(parents=True, exist_ok=True)
    register_tree_program()
    if rule != current:
        write_model(work, rule)

    attempt = {
        "n": len(settings["attempts"]) + 1, "resume": resume and not restart_reason, "command": cmd,
        "environment": env_add, "search_restarted": restart_reason,
        "versions": versions, "code": code, "coordinator_version": cichlid_data.__version__,
        "settings": {"method_options": METHOD, "tree_program": TREE_ENTRY, "iqtree_model": rule,
                     "model_settings": str(work / LIVE_SETTINGS),
                     "seed": seed, "thread_rule": thread_rule, "cores": cores, "tree_order": order,
                     "launcher": "cichlid-orthofinder-of (OrthoFinder's own cmd_order option, set without changing its files)",
                     "launcher_patch": {"smallest_first_fix_in_use": order == "ascending",
                                        "trees_msa_md5": TREES_MSA_MD5},
                     "model": {"value": model, "rule": rule},
                     "iqtree_argument_template": code["iqtree_argument_template"],
                     "species_tree": str(tree) if tree else None,
                     "threads": threads, "analysis_threads": analysis_threads,
                     "snapshot_sizes": list(snapshot_sizes)},
        "host": os.uname().nodename, "started": _now(),
    }
    settings["attempts"].append(attempt)
    settings_file.parent.mkdir(parents=True, exist_ok=True)
    settings_file.write_text(json.dumps(settings, indent=2) + "\n")

    def save():
        settings_file.write_text(json.dumps(settings, indent=2) + "\n")

    watcher = None
    if snapshot_sizes:
        from .snapshot import Watcher
        watcher = Watcher(work, output_dir, snapshot_sizes, tree, state=settings.setdefault("interim", {}), say=say)

    say("$ " + " ".join(cmd))
    t0 = time.time()
    log_file = work / f"attempt_{attempt['n']}.log"
    rc = run_watched(cmd, {**os.environ, **env_add}, log_file, watcher, poll_seconds, on_snapshot=save, say=say)
    attempt.update(finished=_now(), wall_seconds=round(time.time() - t0), exit_code=rc)
    ok = finished_ok(rc, log_file)
    if ok:
        missing = missing_trees(latest_results(output_dir))
        if missing:                                          # OrthoFinder "finished" without some trees
            ok = False
            attempt["missing_trees"] = missing
            say(f"OrthoFinder wrote its results without {len(missing)} gene tree(s) "
                f"(e.g. {', '.join(f'{c} ({n} sequences)' for c, n in missing[:3])}): not accepted")
    if ok:
        attempt["results"] = str(latest_results(output_dir))
        settings["results"] = attempt["results"]
        if watcher:
            watcher.finish(attempt["results"])
    save()
    if not ok:
        miss = attempt.get("missing_trees")
        raise RuntimeError(f"OrthoFinder did not finish cleanly (exit {rc}"
                           + (f"; {len(miss)} clusters without a tree" if miss else "")
                           + f"); see {log_file}; rerun with --resume (finished trees come from the cache)")
    say(f"done in {attempt['wall_seconds']} s; results in {attempt['results']}; settings in {settings_file}")
    return settings_file


def finished_ok(rc, log_file, results=None):
    """Success: exit 0, OrthoFinder's completion message, no traceback, no failed external program, and (when
    the Results folder is given) a tree for every cluster of 3+ sequences. OrthoFinder can exit 0 after a crash
    and after dropping failed trees (docs/pilot.md), so its exit code alone means nothing."""
    text = ANSI.sub("", Path(log_file).read_text(errors="replace"))
    ok = (rc == 0 and ("OrthoFinder run completed" in text or "Done writing files" in text)
          and "Traceback" not in text and "ERROR: external program returned code" not in text)
    if ok and results is not None:
        ok = not missing_trees(results)
    return ok


def missing_trees(results):
    """-> [(cluster, sequences)] for clusters of 3+ sequences without a tree. OrthoFinder carries on and
    reports success when a tree job fails (docs/pilot.md), so a finished run is checked for this."""
    wd = Path(results) / "WorkingDirectory"
    have = {p.stem for p in (wd / "Trees_ids").glob("*.txt")}
    out = []
    for f in sorted((wd / "Sequences_ids").glob("*.fa")):
        data = f.read_bytes()
        n = data.count(b"\n>") + data.startswith(b">")
        if n >= 3 and f.stem not in have:
            out.append((f.stem, n))
    return sorted(out, key=lambda x: -x[1])


def run_watched(cmd, env, log_file, watcher=None, poll_seconds=300, on_snapshot=None, say=print):
    """Run OrthoFinder; while it runs, let the watcher take the interim snapshot when the time comes.
    The snapshot works on a copy and never pauses OrthoFinder. -> OrthoFinder's exit code."""
    with open(log_file, "w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
        while True:
            try:
                return proc.wait(timeout=poll_seconds)
            except subprocess.TimeoutExpired:
                pass
            if watcher is not None:
                try:
                    if watcher.poll() and on_snapshot:
                        on_snapshot()
                except Exception as e:                      # a failed snapshot must never stop the run
                    say(f"interim snapshot failed (the run continues): {e}")


def detach(argv, log):
    """Start this command again in its own session, so it survives the terminal closing."""
    args = [a for a in argv if a != "--detach"]
    with open(log, "a") as fh:
        # unbuffered, so runner.log shows each message as it is written
        p = subprocess.Popen([sys.executable, "-m", "cichlid_orthofinder.cli", *args], stdout=fh,
                             stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True,
                             env={**os.environ, "PYTHONUNBUFFERED": "1"})
    return p.pid
