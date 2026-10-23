"""The gene-tree wrapper: OrthoFinder's tree program, with a cache, IQ-TREE checkpoints and a fixed seed.

OrthoFinder runs one command per gene tree, from its configuration (`~/config_orthofinder_user.json`):

    cichlid-orthofinder-tree INPUT PATH/IDENTIFIER METHODTHREAD   ->   PATH/IDENTIFIER.treefile

For each call this wrapper:

1. Reads the alignment and works out the IQ-TREE settings: the model setting and the threads by the number
   of sequences (default: ModelFinder, `-m MFP`, below 1,000 sequences and the shorter `-m TEST` from 1,000),
   and a fixed seed.
2. Keys the work by a SHA-256 of the settings and the alignment's bytes: a cached tree is used only when
   both are identical.
3. In the key's cache folder:
   - finished before ("done.json" present): the cached tree is returned, IQ-TREE is not run;
   - IQ-TREE finished but the wrapper did not record it (killed in between): recorded now, not rerun;
   - started before (IQ-TREE's checkpoint present): IQ-TREE is rerun with the same command and resumes;
     if that fails, the partial files are moved aside and the tree is built afresh;
   - never started: IQ-TREE is run, writing its checkpoint as it goes.
4. Writes the tree to PATH/IDENTIFIER.treefile, atomically, only once it is complete.

Every decision is appended to `<cache>/events.jsonl`; a running job leaves `running.json` in its folder.
Settings come from the environment (set by `cichlid-orthofinder run`, or by hand):

    CICHLID_TREE_CACHE     the cache folder (required)
    CICHLID_TREE_SEED      IQ-TREE's seed (default 1)
    CICHLID_TREE_THREADS   threads by size, e.g. "1;500:8;2000:16" (default "1": one thread per tree)
    CICHLID_TREE_CORES     the most threads all trees together may use (default: the machine's cores)
    CICHLID_TREE_MODEL     IQ-TREE's -m by size, e.g. "MFP;1000:TEST" (default "MFP" for every tree)
    CICHLID_TREE_SETTINGS  a JSON file, {"model": "<rule>"}, read again at the start of every tree: its rule
                           replaces CICHLID_TREE_MODEL, so the model can be changed while the run goes on
                           (`cichlid-orthofinder set-model`). Trees already started keep the setting they had.
    CICHLID_IQTREE         the IQ-TREE program (default "iqtree3")

The thread budget: before running IQ-TREE, a job reserves its threads from a shared pool (files under
`<cache>/budget/`, guarded by a lock); it waits while the pool is short, and releases its threads when it
ends. Reservations of processes that no longer exist are dropped, so a crash cannot leak cores.

The model is part of the cache key: a tree built with MFP has the same key as before the model rule existed,
and a tree switched to another model gets a new key (its MFP checkpoint, if any, is left in place).
"""

import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

FINISHED = "Total wall-clock time used"     # the last line IQ-TREE writes to its log


def _now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def parse_threads(rule):
    """'1;500:8;2000:16' -> [(0, 1), (500, 8), (2000, 16)]: from N sequences up, use T threads."""
    steps = []
    for i, part in enumerate(p.strip() for p in str(rule).split(";") if p.strip()):
        if ":" in part:
            n, t = part.split(":")
            steps.append((int(n), int(t)))
        elif i == 0:
            steps.append((0, int(part)))
        else:
            raise ValueError(f"bad thread rule part: {part!r}")
    if not steps or steps[0][0] != 0:
        steps.insert(0, (0, 1))
    steps.sort()
    if any(t < 1 for _, t in steps):
        raise ValueError(f"threads must be at least 1: {rule!r}")
    return steps


def threads_for(n, steps):
    t = 1
    for start, threads in steps:
        if n >= start:
            t = threads
    return t


MODEL_NAME = re.compile(r"^[A-Za-z0-9+._{}/,-]+$")     # an IQ-TREE -m value: MFP, TEST, LG+G4, ...


def parse_models(rule):
    """'MFP;1000:TEST' -> [(0, 'MFP'), (1000, 'TEST')]: from N sequences up, use IQ-TREE's -m M."""
    steps = []
    for i, part in enumerate(p.strip() for p in str(rule).split(";") if p.strip()):
        if ":" in part:
            n, m = part.split(":", 1)
            steps.append((int(n), m.strip()))
        elif i == 0:
            steps.append((0, part))
        else:
            raise ValueError(f"bad model rule part: {part!r}")
    if not steps or steps[0][0] != 0:
        steps.insert(0, (0, "MFP"))
    steps.sort()
    bad = [m for _, m in steps if not MODEL_NAME.match(m)]
    if bad or len({n for n, _ in steps}) < len(steps):
        raise ValueError(f"bad model rule: {rule!r}")
    return steps


def model_for(n, steps):
    m = "MFP"
    for start, model in steps:
        if n >= start:
            m = model
    return m


def alignment_shape(data):
    """-> (sequences, columns) of an aligned FASTA given as bytes."""
    n, cols, cur = 0, 0, 0
    for line in data.splitlines():
        if line.startswith(b">"):
            if n == 1:
                cols = cur
            n += 1
            cur = 0
        else:
            cur += len(line.strip())
    if n == 1:
        cols = cur
    return n, cols


class Settings:
    def __init__(self, cache, seed=1, threads="1", iqtree="iqtree3", cores=None, model="MFP", live=None):
        self.cache = Path(cache)
        self.seed = int(seed)
        self.rule = str(threads)
        self.steps = parse_threads(threads)
        self.iqtree = iqtree
        self.cores = int(cores) if cores else (os.cpu_count() or 1)
        self.model_rule = str(model)
        parse_models(self.model_rule)
        self.live = Path(live) if live else None

    def models(self):
        """-> (rule, steps) for a tree starting now: the live settings file's rule, if it has a valid one."""
        if self.live is not None and self.live.is_file():
            try:
                rule = str(json.loads(self.live.read_text())["model"])
                return rule, parse_models(rule)
            except (OSError, ValueError, KeyError, TypeError) as e:
                _event(self.cache, event="settings_unreadable", file=str(self.live), error=str(e),
                       using=self.model_rule)
        return self.model_rule, parse_models(self.model_rule)

    @classmethod
    def from_env(cls, env=None):
        env = os.environ if env is None else env
        if not env.get("CICHLID_TREE_CACHE"):
            raise SystemExit("cichlid-orthofinder-tree: CICHLID_TREE_CACHE is not set (the cache folder)")
        return cls(env["CICHLID_TREE_CACHE"], env.get("CICHLID_TREE_SEED", "1"),
                   env.get("CICHLID_TREE_THREADS", "1"), env.get("CICHLID_IQTREE", "iqtree3"),
                   env.get("CICHLID_TREE_CORES"), env.get("CICHLID_TREE_MODEL", "MFP"),
                   env.get("CICHLID_TREE_SETTINGS"))


def _pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except ProcessLookupError:
        return False
    except (PermissionError, ValueError):
        return True


@contextmanager
def reserved(cache, threads, cores, poll=2.0):
    """Hold `threads` of the shared pool of `cores` while the block runs (waits until they are free)."""
    threads = min(threads, cores)
    pool = Path(cache) / "budget"
    pool.mkdir(parents=True, exist_ok=True)
    me = pool / f"{os.getpid()}-{threading.get_ident()}.res"
    while True:
        with _locked(pool / "lock"):
            used = 0
            for f in pool.glob("*.res"):
                pid = f.name.split("-")[0]
                if not _pid_alive(pid):
                    f.unlink(missing_ok=True)          # a dead process's reservation
                    continue
                try:
                    used += int(f.read_text() or 0)
                except (OSError, ValueError):
                    continue
            if used + threads <= cores:
                me.write_text(str(threads))
                break
        time.sleep(poll)
    try:
        yield threads
    finally:
        me.unlink(missing_ok=True)



def iqtree_version(program, runner=subprocess.run):
    out = runner([program, "--version"], capture_output=True, text=True).stdout
    m = re.search(r"IQ-TREE version (\S+)", out or "")
    if not m:
        raise SystemExit(f"cichlid-orthofinder-tree: could not read the version of {program}")
    return m.group(1)


def iqtree_args(aln, prefix, seed, threads, model="MFP"):
    # -m MFP is IQ-TREE's default (the agreed setting), written out so the record says so.
    return ["-s", str(aln), "--prefix", str(prefix), "-m", model, "-seed", str(seed), "-T", str(threads), "-quiet"]


def cache_key(params, alignment):
    h = hashlib.sha256(json.dumps(params, sort_keys=True).encode() + b"\n")
    h.update(alignment)
    return h.hexdigest()


def _atomic_write(path, data):
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "wb") as fh:
        fh.write(data if isinstance(data, bytes) else data.encode())
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _event(cache, **fields):
    line = json.dumps({"time": _now(), **fields}, sort_keys=True) + "\n"
    with open(Path(cache) / "events.jsonl", "a") as fh:
        fh.write(line)


@contextmanager
def _locked(path):
    with open(path, "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def state(folder):
    """'done', 'finished' (IQ-TREE finished, not yet recorded), 'partial' (a checkpoint) or 'new'."""
    folder = Path(folder)
    if (folder / "done.json").is_file() and (folder / "tree.nwk").is_file():
        return "done"
    log = folder / "run.log"
    if (folder / "run.treefile").is_file() and (folder / "run.iqtree").is_file() and log.is_file() \
            and FINISHED in log.read_text(errors="replace"):
        return "finished"
    if any((folder / f"run.{ext}").exists() for ext in ("ckp.gz", "model.gz", "log")):
        return "partial"
    return "new"


def build(input_fn, out_prefix, settings, runner=subprocess.run, version=None):
    """Make PATH/IDENTIFIER.treefile for one alignment. -> the event recorded (a dict)."""
    data = Path(input_fn).read_bytes()
    n, cols = alignment_shape(data)
    threads = threads_for(n, settings.steps)
    settings.cache.mkdir(parents=True, exist_ok=True)
    _, model_steps = settings.models()
    model = model_for(n, model_steps)
    version = version or iqtree_version(settings.iqtree, runner)
    params = {"program": "IQ-TREE", "version": version, "model": model, "seed": settings.seed, "threads": threads}
    key = cache_key(params, data)
    folder = settings.cache / key[:2] / key
    folder.mkdir(parents=True, exist_ok=True)
    ident = Path(out_prefix).name
    common = dict(og=ident, key=key[:16], sequences=n, columns=cols, threads=threads, model=model)

    with _locked(folder / "lock"):
        before = state(folder)
        seconds = 0.0
        if before in ("partial", "new"):
            if not (folder / "aln.fa").is_file():
                _atomic_write(folder / "aln.fa", data)
                _atomic_write(folder / "params.json", json.dumps(params, sort_keys=True, indent=1) + "\n")
            _atomic_write(folder / "running.json",
                          json.dumps({"pid": os.getpid(), "started": _now(), "start_epoch": time.time(), **common}))
            _event(settings.cache, event="resume" if before == "partial" else "start", **common)
            t0 = time.time()
            ok = _run_iqtree(folder, settings, threads, model, runner)
            if not ok and before == "partial":
                # the checkpoint could not be used: keep it aside and build afresh
                aside = folder / f"failed_{int(time.time())}"
                aside.mkdir()
                for f in folder.glob("run.*"):
                    shutil.move(str(f), str(aside / f.name))
                _event(settings.cache, event="restart_fresh", aside=aside.name, **common)
                ok = _run_iqtree(folder, settings, threads, model, runner)
            seconds = round(time.time() - t0, 1)
            (folder / "running.json").unlink(missing_ok=True)
            if not ok:
                _event(settings.cache, event="failed", seconds=seconds, **common)
                raise SystemExit(f"cichlid-orthofinder-tree: IQ-TREE failed for {ident}; see {folder / 'run.log'}")
        if before != "done":
            _atomic_write(folder / "tree.nwk", (folder / "run.treefile").read_bytes())
            _atomic_write(folder / "done.json", json.dumps(
                {"finished": _now(), "params": params, "sequences": n, "columns": cols,
                 "first_og": ident, "how": {"new": "built", "partial": "resumed", "finished": "recorded"}[before]},
                sort_keys=True, indent=1) + "\n")
        tree = (folder / "tree.nwk").read_bytes()

    _atomic_write(Path(f"{out_prefix}.treefile"), tree)
    how = {"done": "reused", "finished": "recorded", "partial": "resumed", "new": "built"}[before]
    ev = dict(event="done", how=how, seconds=seconds, **common)
    _event(settings.cache, **ev)
    return ev


def _run_iqtree(folder, settings, threads, model, runner):
    cmd = [settings.iqtree, *iqtree_args(folder / "aln.fa", folder / "run", settings.seed, threads, model)]
    with reserved(settings.cache, threads, settings.cores):
        runner(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # judged by IQ-TREE's outputs, not its exit code: its log's last line is written only when it is done
    return state(folder) == "finished"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        print("usage: cichlid-orthofinder-tree INPUT PATH/IDENTIFIER [METHODTHREAD]", file=sys.stderr)
        return 2
    build(argv[0], argv[1], Settings.from_env())
    return 0


if __name__ == "__main__":
    sys.exit(main())
