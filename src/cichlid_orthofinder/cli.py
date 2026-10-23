"""The cichlid-orthofinder command."""

import argparse
import os
import sys

from . import __version__


# subcommands whose options belong to their own module: everything after the name goes to its main()
PASSTHROUGH = {
    "convert": ("convert", "OrthoFinder's output -> assign/, families/, run.json and native/summary/"),
    "native": ("native", "the curated copy of OrthoFinder's result folder, for native/"),
    "verify": ("verify", "the gates before a restart, each PASS or FAIL (saved to gates.txt)"),
}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] in PASSTHROUGH:
        import importlib
        module = importlib.import_module(f"cichlid_orthofinder.{PASSTHROUGH[argv[0]][0]}")
        return module.main(argv[1:])
    p = argparse.ArgumentParser(prog="cichlid-orthofinder",
                                description="OrthoFinder families and assignments for the cichlid project.")
    p.add_argument("--version", action="version", version=f"cichlid-orthofinder {__version__}")
    sub = p.add_subparsers(dest="command", required=True)
    for name, (_, help_text) in PASSTHROUGH.items():          # listed in --help; handled above
        sub.add_parser(name, help=help_text + f" (cichlid-orthofinder {name} --help for its options)")

    r = sub.add_parser("run", help="check the inputs and run (or resume) OrthoFinder on the release's species")
    r.add_argument("--work", required=True,
                   help="scratch folder for the input links, OrthoFinder's output and the tree cache (not a data folder)")
    r.add_argument("--release", default="v1")
    r.add_argument("--run", default="our_species_tree", choices=("our_species_tree", "own_species_tree"),
                   help="our_species_tree uses the catalog's tree (-s); own_species_tree lets OrthoFinder infer one")
    r.add_argument("--seed", type=int, default=1, help="IQ-TREE's seed (default 1)")
    r.add_argument("--thread-rule", default="1;600:2;1500:8;3000:16",
                   help="threads per tree by size: '1;600:2;1500:8;3000:16' = 1 thread, 2 from 600 sequences, 8 from 1,500, "
                        "16 from 3,000 (default; chosen by simulating the run order, pilot_runs/scripts/threshold_search.py)")
    r.add_argument("--order", default="ascending", choices=("ascending", "descending"),
                   help="smallest first (default: interim results come early; OrthoFinder's own order, with a one-line in-memory fix that leaves outputs identical, docs/pilot.md) or largest first")
    r.add_argument("--model", default="full",
                   help="IQ-TREE's -m by size: 'full' (MFP for all, the sheet's setting; default), 'hybrid' (MFP below "
                        "1,000 sequences, TEST from 1,000) or a rule such as 'MFP;1000:TEST'")
    r.add_argument("--from-search", help="start from this finished search (a WorkingDirectory with the Blast*.txt.gz "
                                         "files), used only if every file is complete and it is a search of exactly "
                                         "the release's proteomes (otherwise the run searches again); it is copied, the original "
                                         "is not touched")
    r.add_argument("--cores", type=int, help="most threads all trees together may use (default: all cores)")
    r.add_argument("--threads", type=int, help="OrthoFinder -t (parallel jobs; default: --cores)")
    r.add_argument("--analysis-threads", type=int, help="OrthoFinder -a")
    r.add_argument("--snapshot-sizes", default="300,1000",
                   help="interim results once every cluster of up to each of these sizes has its tree "
                        "(default '300,1000'; '' = none); the final results come when the run ends")
    r.add_argument("--resume", action="store_true",
                   help="restart an interrupted run from its search; finished trees come from the cache")
    r.add_argument("--detach", action="store_true", help="run in the background, surviving logout")
    r.add_argument("--dry-run", action="store_true", help="check the inputs and print the command only")

    sm = sub.add_parser("set-model", help="change the model setting of a running run: trees that start from now "
                                          "use it; trees already running or finished keep theirs")
    sm.add_argument("work", help="the run's work folder (<work>/<release>/<run>)")
    sm.add_argument("model", nargs="?", help="'hybrid', 'full' or a rule such as 'MFP;1000:TEST' "
                                              "(leave out to show the current setting)")

    s = sub.add_parser("status", help="progress of a running or finished OrthoFinder run (reads only)")
    s.add_argument("results", help="the run's work folder (<work>/<release>/<run>), its output folder, "
                                   "or a Results_* folder")
    s.add_argument("--cache", default=os.environ.get("CICHLID_TREE_CACHE"),
                   help="the tree wrapper's cache folder (default: <work>/<release>/<run>/tree_cache if found)")
    s.add_argument("--threads", default=os.environ.get("CICHLID_TREE_THREADS"),
                   help="the wrapper's thread rule (default: the run's, from run_settings.json)")

    n = sub.add_parser("snapshot", help="interim results from the finished trees of a running run (a copy; "
                                        "the run is untouched)")
    n.add_argument("results", help="the run's work folder (<work>/<release>/<run>) or a Results_* folder")
    n.add_argument("dest", help="a new folder for the snapshot")
    n.add_argument("--species-tree", help="the -s species tree (default: the run's, from run_settings.json)")
    n.add_argument("--threads", type=int, default=8)
    m = sub.add_parser("map", help="snapshot orthogroup -> final orthogroup, once the final run is done")
    m.add_argument("snapshot", help="the snapshot's Results_snapshot_* folder")
    m.add_argument("final", help="the final Results_* folder")
    m.add_argument("out", help="the TSV to write")

    args = p.parse_args(argv)
    if args.command == "snapshot":
        from pathlib import Path
        import json
        from .run import latest_results
        from .snapshot import take
        path, tree = Path(args.results).resolve(), args.species_tree
        settings = path / "run_settings.json"
        if settings.is_file():
            last = json.loads(settings.read_text())["attempts"][-1]
            tree = tree or last["settings"]["species_tree"]
            path = path / "output"
        if not path.name.startswith("Results_"):
            path = latest_results(path)
        take(path, args.dest, tree, args.threads)
        return 0
    if args.command == "set-model":
        from pathlib import Path
        from .run import MODEL_SHORTCUTS, read_model, write_model
        work = Path(args.work).expanduser().resolve()
        if not (work / "run_settings.json").is_file():
            print(f"error: not a run's work folder (no run_settings.json): {work}", file=sys.stderr)
            return 1
        names = {v: k for k, v in MODEL_SHORTCUTS.items()}
        if args.model:
            try:
                write_model(work, args.model)
            except ValueError as e:
                print(f"error: {e}", file=sys.stderr)
                return 1
        rule = read_model(work)
        print(f"model setting: {rule}" + (f" ({names[rule]})" if rule in names else "")
              + ("; trees that start from now use it" if args.model else ""))
        return 0
    if args.command == "map":
        from .snapshot import mapping
        s = mapping(args.snapshot, args.final, args.out)
        print(f"{s['matched']} of {s['orthogroups']} snapshot orthogroups found in the final run "
              f"({s['same_tree']} with the same resolved tree); written to {args.out}")
        return 0
    if args.command == "status":
        from pathlib import Path
        import json
        from .run import latest_results
        from .status import report
        path = Path(args.results).resolve()
        from datetime import datetime
        from .run import read_model
        cache, rule, cores, model, started = args.cache, args.threads, None, "MFP", None
        settings = path / "run_settings.json"
        if settings.is_file():                       # a run's work folder
            s = json.loads(settings.read_text())
            last = s["attempts"][-1]
            cache = cache or last["environment"]["CICHLID_TREE_CACHE"]
            rule = rule or last["settings"]["thread_rule"]
            cores = last["settings"].get("cores")
            model = read_model(path) or last["settings"].get("iqtree_model", "MFP")
            started = datetime.fromisoformat(s["attempts"][0]["started"]).replace(tzinfo=None)
            path = path / "output"
        if not path.name.startswith("Results_"):
            path = latest_results(path)
        report(path, cache=cache, threads_rule=rule or "1", cores=cores, model_rule=model, started=started)
        return 0
    if args.command == "run":
        from .run import detach, run
        if args.detach:
            from pathlib import Path
            log = Path(args.work).expanduser().resolve() / args.release / args.run
            log.mkdir(parents=True, exist_ok=True)
            pid = detach(argv, log / "runner.log")
            print(f"running in the background (pid {pid}); log: {log / 'runner.log'}")
            print(f"progress: cichlid-orthofinder status {log}")
            return 0
        try:
            run(args.work, run_name=args.run, release=args.release, threads=args.threads,
                analysis_threads=args.analysis_threads, seed=args.seed, thread_rule=args.thread_rule,
                cores=args.cores, order=args.order, model=args.model, from_search=args.from_search, snapshot_sizes=[int(s) for s in args.snapshot_sizes.split(',') if s.strip()],
                resume=args.resume, dry_run=args.dry_run)
        except (RuntimeError, FileNotFoundError, ValueError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
