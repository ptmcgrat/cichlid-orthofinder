"""Progress of a running (or finished) OrthoFinder run, read from its result folder and, when given, the
tree wrapper's cache (treecache.py). Reads only.

Stages, as OrthoFinder 3.1.5 writes them (docs/pilot.md):
    search      WorkingDirectory/Blast<i>_<j>.txt.gz appear, one per ordered species pair
    clustering  clusters_OrthoFinder_I*.txt
    trees       one Sequences_ids/<cluster>.fa per preliminary cluster, then Trees_ids/<cluster>.txt as
                each tree finishes (clusters of 3+ sequences get one)
    final       Resolved_Gene_Trees/, Gene_Duplication_Events/ and the other result folders
    done        "OrthoFinder run completed" in Log.txt

Tree progress is measured in expected work, not in trees. The cost model is the plan's (docs/pilot.md,
pilot_runs/scripts/schedule_sim.py): one IQ-TREE run on n sequences takes about C * n^B seconds on one
thread up to BIG sequences; above that, GIANT_ANCHOR seconds at ANCHOR_N sequences (the pilot's largest real
tree) times (n / ANCHOR_N)^GIANT_K, halved when the model setting for that size is TEST (assumed, not
measured). T threads make a tree SPEED[T] times faster at CPU_COST[T] times the CPU. Once trees of this run
have finished, the model is scaled by how long they actually took (the median of actual / expected, from
trees of CALIBRATE_FROM sequences or more once there are enough: the smallest are mostly fixed overhead).

The time to go is the larger of two numbers: the remaining work at the pace measured so far, and the
longest time any single remaining tree still needs on its own threads. The dates of the three stages
(every cluster up to 300 sequences, up to 1,000, all) come from simulating the remaining trees in
OrthoFinder's order (smallest first) on the run's cores with the thread rule, the running trees first.
"""

import json
import os
import re
import statistics
import subprocess
import time
from pathlib import Path

MIN_TREE = 3            # clusters with fewer sequences get no tree (pilot: every cluster of 3+ had one)
C, B = 0.148 * 0.84, 2.23   # seconds on one thread for n sequences: C * n**B (pilot fit, scaled to its total)
SPEED = {1: 1.0, 2: 2.0, 4: 2.5, 8: 3.3, 16: 4.4}       # measured on the pilot's largest tree
CPU_COST = {1: 1.0, 2: 1.0, 4: 1.4, 8: 1.9, 16: 2.4}
BIG, ANCHOR_N, GIANT_ANCHOR, GIANT_K = 1000, 576, 12.5 * 3600, 1.5   # the plan's model above 1,000 sequences
CALIBRATE_FROM, CALIBRATE_MIN = 20, 30
STAGES = (300, 1000)


def _nearest(table, t):
    return table[max(k for k in table if k <= t)] if t >= 1 else 1.0


def single_seconds(n, model="MFP"):
    """Expected time of one tree of n sequences on one thread (before calibration)."""
    if n <= BIG:
        return C * n ** B
    return GIANT_ANCHOR * (n / ANCHOR_N) ** GIANT_K * (0.5 if model == "TEST" else 1.0)


def expected_seconds(n, threads=1, factor=1.0, model="MFP"):
    """Expected wall time of one tree of n sequences on `threads` threads."""
    return factor * single_seconds(n, model) / _nearest(SPEED, threads)


def core_seconds(n, threads=1, factor=1.0, model="MFP"):
    """Expected CPU (core) time of one tree: the unit of work."""
    return factor * single_seconds(n, model) * _nearest(CPU_COST, threads)


def simulate(queue, running_left, cores):
    """Finish times (seconds from now) of trees started in `queue` order [(cluster, n, seconds on its
    threads, threads)] once `running_left` [(cluster, n, seconds left, threads)] are under way, each tree
    waiting until its threads are free (the wrapper's shared budget). -> {cluster: (n, finish)}."""
    import heapq
    t, free, heap, finish = 0.0, cores, [], {}
    for og, n, left, th in running_left:
        th = min(th, cores)
        heapq.heappush(heap, (left, th, og, n))
        free -= th
    for og, n, dur, th in queue:
        th = min(th, cores)
        while free < th and heap:
            end, c, o, m = heapq.heappop(heap)
            t = max(t, end); free += c; finish[o] = (m, end)
        heapq.heappush(heap, (t + dur, th, og, n))
        free -= th
    while heap:
        end, c, o, m = heapq.heappop(heap)
        finish[o] = (m, end)
    return finish


def _count_fasta(path):
    data = Path(path).read_bytes()
    return data.count(b"\n>") + data.startswith(b">")


def cluster_sizes(working):
    """-> {cluster: number of sequences} for every preliminary cluster with a sequence file."""
    folder = Path(working) / "Sequences_ids"
    if not folder.is_dir():
        return {}
    return {p.stem: _count_fasta(p) for p in folder.glob("*.fa")}


def running_iqtree(working, runner=subprocess.run):
    """-> [(cluster, elapsed seconds)] for IQ-TREE processes running directly on this run's alignments
    (runs without the wrapper)."""
    try:
        out = runner(["ps", "-eo", "etimes=,args="], capture_output=True, text=True).stdout
    except OSError:
        return []
    mine = re.compile(re.escape(str(Path(working) / "Alignments_ids")) + r"/(OG\d+)\.fa\b")
    jobs = []
    for line in out.splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2 and parts[1].startswith("iqtree3 "):
            m = mine.search(parts[1])
            if m:
                jobs.append((m.group(1), int(parts[0])))
    return jobs


LEFT = re.compile(r"Iteration (\d+) / LogL: \S+ / Time: \S+ \((\d+)h:(\d+)m:(\d+)s left\)")


def job_state(log):
    """-> (what, seconds_left) for one running IQ-TREE job, from its log. seconds_left is IQ-TREE's own
    estimate during the tree search (it restarts whenever a better tree is found), else None."""
    try:
        text = Path(log).read_text(errors="replace")
    except OSError:
        return "?", None
    if "Best-fit model" not in text:
        return "choosing a model", None
    found = LEFT.findall(text)
    if found:
        it, h, m, s = (int(x) for x in found[-1])
        return f"tree search, iteration {it}", h * 3600 + m * 60 + s
    if "TREE SEARCH COMPLETED" in text:
        return "finishing", 0
    return "starting the tree search", None


def models_tested(model_gz, runner=subprocess.run):
    """Models ModelFinder has finished so far (one line each in IQ-TREE's <prefix>.model.gz)."""
    if not Path(model_gz).is_file():
        return 0
    try:
        out = runner(["zcat", str(model_gz)], capture_output=True, text=True).stdout
    except OSError:
        return 0
    return sum(1 for l in out.splitlines() if ":" in l and not l.startswith("---"))


def running_diamond(working, runner=subprocess.run):
    """-> [pair file] for the DIAMOND searches running now for this run."""
    try:
        out = runner(["ps", "-eo", "args="], capture_output=True, text=True).stdout
    except OSError:
        return []
    mine = re.compile(re.escape(str(Path(working))) + r"/(Blast\d+_\d+\.txt)")
    return [m.group(1) for line in out.splitlines()
            if line.startswith("diamond blastp") and (m := mine.search(line))]


def stage(results, runner=subprocess.run):
    """-> (stage, detail) for a Results_* folder."""
    results = Path(results)
    w = results / "WorkingDirectory"
    log = results / "Log.txt"
    if log.is_file() and "OrthoFinder run completed" in log.read_text(errors="replace"):
        return "done", ""
    if (results / "Resolved_Gene_Trees").is_dir():
        return "final", "writing the result folders"
    if (w / "Sequences_ids").is_dir():
        return "trees", ""
    if list(w.glob("clusters_OrthoFinder_I*.txt")):
        return "clustering", ""
    n_species = len(list(w.glob("Species*.fa")))
    blast = list(w.glob("Blast*_*.txt.gz"))
    if n_species:
        # DIAMOND writes each pair's file as it goes: a pair still being searched has a file too
        busy = len(running_diamond(w, runner))
        return "search", f"{max(len(blast) - busy, 0)} of {n_species ** 2} species pairs searched, {busy} running"
    return "starting", ""


def _alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def cache_state(cache, now=None):
    """-> (events summary, running jobs, calibration samples) from the wrapper's cache."""
    cache = Path(cache)
    hows, failed, restarted, samples = {}, 0, 0, []
    ev_file = cache / "events.jsonl"
    if ev_file.is_file():
        for line in ev_file.read_text().splitlines():
            try:
                e = json.loads(line)
            except json.JSONDecodeError:
                continue
            if e.get("event") == "done":
                hows[e["how"]] = hows.get(e["how"], 0) + 1
                if e["how"] == "built" and e.get("seconds", 0) > 2:
                    samples.append((e["sequences"], e["seconds"] / expected_seconds(
                        e["sequences"], e["threads"], model=e.get("model", "MFP"))))
            elif e.get("event") == "failed":
                failed += 1
            elif e.get("event") == "restart_fresh":
                restarted += 1
    now = now or time.time()
    running = []
    for rj in cache.glob("*/*/running.json"):
        try:
            r = json.loads(rj.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not _alive(r.get("pid")):
            continue
        folder = rj.parent
        what, left = job_state(folder / "run.log")
        if what == "choosing a model":
            what = f"choosing a model ({models_tested(folder / 'run.model.gz')} tested)"
        running.append(dict(og=r["og"], sequences=r["sequences"], threads=r["threads"],
                            elapsed=now - r["start_epoch"], what=what, iqtree_left=left))
    return dict(hows=hows, failed=failed, restarted=restarted), running, samples


def progress(results, cache=None, threads_rule="1", now=None, jobs=None, cores=None, model_rule="MFP"):
    """-> dict: stage, counts, work done, time to go, stage dates, running jobs (with the wrapper's cache
    if given)."""
    from .treecache import parse_models, parse_threads, model_for, threads_for
    results = Path(results)
    w = results / "WorkingDirectory"
    st, detail = stage(results)
    out = {"results": str(results), "stage": st, "detail": detail}
    if st not in ("trees", "final", "done"):
        return out
    now = now or time.time()
    steps = parse_threads(threads_rule)
    msteps = parse_models(model_rule)
    cores = int(cores or os.cpu_count() or 1)

    sizes = {k: v for k, v in cluster_sizes(w).items() if v >= MIN_TREE}
    finished = {p.stem for p in (w / "Trees_ids").glob("*.txt")} if (w / "Trees_ids").is_dir() else set()

    summary, running, samples = ({}, [], []) if cache is None else cache_state(cache, now)
    if cache is None:
        if jobs is None:
            jobs = running_iqtree(w)
        for og, sec in jobs:
            what, left = job_state(w / "Alignments_ids" / f"{og}.log")
            running.append(dict(og=og, sequences=sizes.get(og), threads=1, elapsed=sec, what=what, iqtree_left=left))
    big = [r for n, r in samples if n >= CALIBRATE_FROM]
    use = big if len(big) >= CALIBRATE_MIN else [r for _, r in samples]
    factor = statistics.median(use) if len(use) >= 5 else 1.0

    def th(og):
        return threads_for(sizes[og], steps)

    def mod(og):
        return model_for(sizes[og], msteps)

    def work(og):
        return core_seconds(sizes[og], th(og), factor, mod(og))

    total = sum(work(og) for og in sizes)
    done = sum(work(og) for og in finished if og in sizes)
    partial, floor = 0.0, 0.0
    run_ogs = set()
    for j in running:
        og, n = j["og"], j["sequences"]
        if og not in sizes or og in finished:
            continue
        run_ogs.add(og)
        exp = expected_seconds(n, j["threads"], factor, model_for(n, msteps))
        j["expected_left"] = max(exp - j["elapsed"], 0.0)
        if j["iqtree_left"]:
            j["expected_left"] = max(j["expected_left"], j["iqtree_left"])
        partial += work(og) * min(j["elapsed"] / exp, 0.95) if exp else 0.0
        floor = max(floor, j["expected_left"])
    waiting = sorted((og for og in sizes if og not in finished and og not in run_ogs), key=lambda o: (sizes[o], o))
    for og in waiting:
        floor = max(floor, expected_seconds(sizes[og], th(og), factor, mod(og)))

    # the three stages: when the last cluster up to each size is expected to finish
    finish = simulate([(og, sizes[og], expected_seconds(sizes[og], th(og), factor, mod(og)), th(og)) for og in waiting],
                      [(j["og"], sizes[j["og"]], j["expected_left"], j["threads"]) for j in running if j["og"] in run_ogs],
                      cores)
    stages = {}
    for limit in (*STAGES, None):
        left = [f for n, f in finish.values() if limit is None or n <= limit]
        stages["all" if limit is None else limit] = max(left) / 3600 if left else 0.0

    out.update(clusters=len(sizes), trees=len(finished & set(sizes)), work_done=(done + partial) / total if total else 0.0,
               running=sorted(running, key=lambda j: -j["elapsed"]), cache=summary, calibrated=len(use),
               factor=factor, stage_hours=stages, cores=cores, model_rule=model_rule,
               stage_counts={lim: (sum(1 for o in sizes if sizes[o] <= lim and o in finished),
                                   sum(1 for o in sizes if sizes[o] <= lim)) for lim in STAGES})
    seq_dir = w / "Sequences_ids"
    started = min((p.stat().st_mtime for p in seq_dir.glob("*.fa")), default=None)
    if st == "trees" and started:
        elapsed = max(now - started, 1.0)
        pace = (done + partial) / elapsed                    # core-seconds of work per second (effective cores)
        remaining = total - done - partial
        by_pace = remaining / pace if pace > 0 else None
        out.update(elapsed_hours=elapsed / 3600, effective_cores=pace,
                   pace_hours=None if by_pace is None else by_pace / 3600, floor_hours=floor / 3600,
                   eta_hours=None if by_pace is None else max(by_pace, floor) / 3600)
    return out


def bar(fraction, width=30):
    full = int(round(fraction * width))
    return "[" + "#" * full + "." * (width - full) + "]"


def _h(hours):
    return f"{hours * 60:.0f} min" if hours < 1 else (f"{hours:.1f} h" if hours < 48 else f"{hours / 24:.1f} days")


def report(results, cache=None, threads_rule="1", say=print, cores=None, model_rule="MFP", started=None):
    p = progress(results, cache, threads_rule, cores=cores, model_rule=model_rule)
    say(f"results: {p['results']}")
    say(f"stage:   {p['stage']}" + (f" ({p['detail']})" if p["detail"] else ""))
    if "clusters" not in p:
        return p
    say(f"work:    {bar(p['work_done'])} {100 * p['work_done']:.1f}%  "
        f"(trees {p['trees']:,} of {p['clusters']:,}; {100 * p['trees'] / max(p['clusters'], 1):.1f}% by count)")
    if p["cache"]:
        c = p["cache"]["hows"]
        say(f"cache:   built {c.get('built', 0):,}  resumed {c.get('resumed', 0):,}  reused {c.get('reused', 0):,}  "
            f"recorded {c.get('recorded', 0):,}  failed {p['cache']['failed']}  restarted fresh {p['cache']['restarted']}")
    if "eta_hours" in p and p["eta_hours"] is not None:
        say(f"elapsed: {_h(p['elapsed_hours'])} in the tree stage, at {p['effective_cores']:.1f} cores' worth of work")
        final = p.get("stage_hours", {}).get("all")
        if final is not None:
            say(f"to go:   about {_h(final)} (simulated in run order on {p['cores']} cores); at the pace so far "
                f"{_h(p['eta_hours'])}, pessimistic while the smallest trees run (mostly fixed overhead)")
        else:
            say(f"to go:   about {_h(p['eta_hours'])} (remaining work at this pace {_h(p['pace_hours'])}; "
                f"slowest remaining tree alone {_h(p['floor_hours'])})")
        say(f"model:   the plan's tree cost, " + (f"calibrated on {p['calibrated']} trees of this run (x{p['factor']:.2f})"
                                                  if p["calibrated"] >= 5 else "not yet calibrated on this run (pilot values)")
            + f"; model setting {p['model_rule']}")
    if p.get("stage_hours"):
        from datetime import datetime, timedelta
        now = datetime.now()
        names = {300: "early results 1 (up to 300 genes)", 1000: "early results 2 (up to 1,000 genes)", "all": "final results"}
        for k, hours in p["stage_hours"].items():
            when = now + timedelta(hours=hours)
            day = f", day {(when - started).total_seconds() / 86400:.1f}" if started else ""
            got = p["stage_counts"].get(k) if k != "all" else (p["trees"], p["clusters"])
            say(f"stage:   {names[k]:36s} ~{when:%a %b %d %H:%M}{day}  ({got[0]:,} of {got[1]:,} trees)"
                if hours > 0 else f"stage:   {names[k]:36s} ready  ({got[0]:,} of {got[1]:,} trees)")
    if p["running"]:
        say(f"running: {len(p['running'])} tree(s); longest first:")
        for j in p["running"][:8]:
            left = f"; ~{_h(j['expected_left'] / 3600)} left" if j.get("expected_left") is not None else ""
            say(f"  {j['og']}  {j['sequences'] if j['sequences'] is not None else '?'} seq  {j['threads']} thr  "
                f"{_h(j['elapsed'] / 3600)}  {j['what']}{left}")
    return p
