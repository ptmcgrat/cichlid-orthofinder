import json
import os
import subprocess
import tempfile
import time
import unittest
import unittest.mock
from pathlib import Path

from cichlid_orthofinder import status


def fasta(path, n):
    path.write_text("".join(f">0_{i}\nMKV\n" for i in range(n)))


class FakeRun:
    """A Results_* folder shaped like OrthoFinder 3.1.5's (docs/pilot.md)."""

    def __init__(self, root):
        self.r = Path(root) / "Results_Oct20"
        self.w = self.r / "WorkingDirectory"
        self.w.mkdir(parents=True)
        (self.r / "Log.txt").write_text("Started OrthoFinder version 3.1.5\n")

    def species(self, n):
        for i in range(n):
            (self.w / f"Species{i}.fa").write_text(">0_0\nM\n")

    def blast(self, pairs):
        for i, j in pairs:
            (self.w / f"Blast{i}_{j}.txt.gz").write_bytes(b"")

    def clusters(self, sizes, mtime=None):
        (self.w / "clusters_OrthoFinder_I1.2.txt").write_text("")
        (self.w / "Sequences_ids").mkdir()
        (self.w / "Trees_ids").mkdir()
        for og, n in sizes.items():
            p = self.w / "Sequences_ids" / f"{og}.fa"
            fasta(p, n)
            if mtime:
                os.utime(p, (mtime, mtime))

    def tree(self, og):
        (self.w / "Trees_ids" / f"{og}.txt").write_text("(0_0,0_1,0_2);")


class TestStage(unittest.TestCase):
    def test_stages_in_order(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            self.assertEqual(status.stage(f.r)[0], "starting")
            f.species(3)
            f.blast([(0, 0), (0, 1)])
            # one of the two pair files is still being written by a running DIAMOND search
            ps = (f"diamond blastp --ignore-warnings -d {f.w}/diamondDBSpecies1 -q {f.w}/Species0.fa "
                  f"-o {f.w}/Blast0_1.txt --matrix BLOSUM62\n"
                  "diamond blastp -d /other/run/WorkingDirectory/x -o /other/run/WorkingDirectory/Blast0_0.txt\n")

            def fake(cmd, **kw):
                return subprocess.CompletedProcess(cmd, 0, stdout=ps, stderr="")

            self.assertEqual(status.stage(f.r, fake), ("search", "1 of 9 species pairs searched, 1 running"))
            f.clusters({"OG0000000": 5})
            self.assertEqual(status.stage(f.r)[0], "trees")
            (f.r / "Resolved_Gene_Trees").mkdir()
            self.assertEqual(status.stage(f.r)[0], "final")
            (f.r / "Log.txt").write_text("...\n2026-10-20 14:43:09 : OrthoFinder run completed\n")
            self.assertEqual(status.stage(f.r)[0], "done")


class TestProgress(unittest.TestCase):
    def test_work_weighted_by_the_cost_model(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            f.species(2)
            # clusters below 3 sequences get no tree
            f.clusters({"OG0000000": 10, "OG0000001": 5, "OG0000002": 5, "OG0000003": 3, "OG0000004": 2},
                       mtime=1000)
            f.tree("OG0000000")
            p = status.progress(f.r, now=1000 + 3600, jobs=[])
            self.assertEqual((p["clusters"], p["trees"]), (4, 1))
            w = {n: status.core_seconds(n) for n in (10, 5, 3)}
            self.assertAlmostEqual(p["work_done"], w[10] / (w[10] + 2 * w[5] + w[3]))
            # pace: w[10] core-seconds in one hour; the rest at that pace
            self.assertAlmostEqual(p["pace_hours"], (2 * w[5] + w[3]) / w[10])

    def test_time_to_go_is_at_least_the_slowest_remaining_tree(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            f.clusters({**{f"OG{i:07d}": 5 for i in range(1, 200)}, "OG0000000": 2000}, mtime=1000)
            for i in range(1, 200):
                f.tree(f"OG{i:07d}")
            p = status.progress(f.r, now=1000 + 60, jobs=[])
            alone = status.expected_seconds(2000) / 3600
            self.assertAlmostEqual(p["floor_hours"], alone)
            self.assertAlmostEqual(p["eta_hours"], max(p["pace_hours"], alone))
            # with 16 threads for clusters of 2,000+, the giant alone takes 4.4 times less
            p16 = status.progress(f.r, threads_rule="1;2000:16", now=1000 + 60, jobs=[])
            self.assertAlmostEqual(p16["floor_hours"], alone / 4.4)

    def test_running_jobs_without_cache(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            f.clusters({"OG0000000": 10, "OG0000001": 5})
            p = status.progress(f.r, jobs=[("OG0000001", 60), ("OG0000000", 7200)])
            self.assertEqual([(j["og"], j["sequences"], j["elapsed"]) for j in p["running"]],
                             [("OG0000000", 10, 7200), ("OG0000001", 5, 60)])


class TestCache(unittest.TestCase):
    def make_cache(self, d, events, running=()):
        cache = Path(d) / "cache"
        cache.mkdir()
        (cache / "events.jsonl").write_text("".join(json.dumps(e) + "\n" for e in events))
        for i, (og, n, pid, start, log) in enumerate(running):
            folder = cache / "ab" / f"key{i}"
            folder.mkdir(parents=True)
            (folder / "running.json").write_text(json.dumps(
                dict(pid=pid, og=og, sequences=n, threads=1, start_epoch=start)))
            (folder / "run.log").write_text(log)
        return cache

    def test_summary_calibration_and_running(self):
        with tempfile.TemporaryDirectory() as d:
            # five trees of 100 sequences that took twice the model's time -> factor 2
            t = status.expected_seconds(100)
            events = [dict(event="done", how="built", sequences=100, threads=1, seconds=2 * t)] * 5
            events += [dict(event="done", how="reused", sequences=50, threads=1, seconds=0),
                       dict(event="failed"), dict(event="restart_fresh")]
            cache = self.make_cache(d, events, running=[
                ("OG0000000", 10, os.getpid(), 1000.0, TestJobState.SEARCH),        # alive
                ("OG0000009", 10, 999999999, 1000.0, "")])                          # dead pid: ignored
            summary, running, samples = status.cache_state(cache, now=1600.0)
            self.assertEqual(summary["hows"], {"built": 5, "reused": 1})
            self.assertEqual((summary["failed"], summary["restarted"]), (1, 1))
            self.assertEqual([round(r, 6) for _, r in samples], [2.0] * 5)
            self.assertEqual([(j["og"], j["elapsed"]) for j in running], [("OG0000000", 600.0)])
            self.assertEqual(running[0]["what"], "tree search, iteration 380")

    def test_report_with_cache(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            f.clusters({"OG0000000": 10, "OG0000001": 5}, mtime=1000)
            f.tree("OG0000001")
            cache = self.make_cache(d, [dict(event="done", how="built", sequences=5, threads=1, seconds=3)],
                                    running=[("OG0000000", 10, os.getpid(), time.time() - 30, TestJobState.MODEL)])
            lines = []
            status.report(f.r, cache=cache, say=lines.append)
            text = "\n".join(lines)
            self.assertIn("cache:   built 1", text)
            self.assertIn("OG0000000  10 seq  1 thr", text)
            self.assertIn("choosing a model (0 tested)", text)
            self.assertIn("work:    [", text)

    def test_bar(self):
        self.assertEqual(status.bar(0.5, 10), "[#####.....]")
        self.assertEqual(status.bar(0.0, 4), "[....]")


class TestJobState(unittest.TestCase):
    # Lines as IQ-TREE 3.1.3 writes them (pilot logs).
    MODEL = ("Alignment has 576 sequences with 525 columns, 525 distinct patterns\n"
             "ModelFinder will test up to 1232 protein models (sample size: 525 epsilon: 0.100) ...\n"
             " No. Model         -LnL         df  AIC          AICc         BIC\n")
    SEARCH = (MODEL + "Best-fit model: JTT+I+R7 chosen according to BIC\n"
              "Iteration 370 / LogL: -83492.816 / Time: 7h:45m:57s (0h:44m:30s left)\n"
              "BETTER TREE FOUND at iteration 375: -83490.1\n"
              "Iteration 380 / LogL: -83496.648 / Time: 7h:50m:19s (0h:39m:10s left)\n")
    DONE = SEARCH + "TREE SEARCH COMPLETED AFTER 476 ITERATIONS / Time: 8h:30m:1s\n"

    def state(self, text):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "OG.log"
            p.write_text(text)
            return status.job_state(p)

    def test_choosing_a_model(self):
        self.assertEqual(self.state(self.MODEL), ("choosing a model", None))

    def test_tree_search_uses_last_estimate(self):
        self.assertEqual(self.state(self.SEARCH), ("tree search, iteration 380", 39 * 60 + 10))

    def test_missing_log(self):
        self.assertEqual(status.job_state("/nonexistent/OG.log"), ("?", None))

    def test_running_job_uses_iqtree_estimate_when_larger(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            f.clusters({"OG0000000": 10, "OG0000001": 5}, mtime=1000)
            f.tree("OG0000001")
            (f.w / "Alignments_ids").mkdir()
            (f.w / "Alignments_ids" / "OG0000000.log").write_text(self.SEARCH)
            lines = []
            with unittest.mock.patch.object(status, "running_iqtree", return_value=[("OG0000000", 3600)]):
                status.report(f.r, say=lines.append)
            text = "\n".join(lines)
            # the model expects seconds for 10 sequences; IQ-TREE says 39 min 10 s: the larger wins
            self.assertIn("tree search, iteration 380; ~39 min left", text)

    def test_before_trees_has_no_counts(self):
        with tempfile.TemporaryDirectory() as d:
            f = FakeRun(d)
            f.species(2)
            self.assertNotIn("clusters", status.progress(f.r, jobs=[]))


class TestRunningIqtree(unittest.TestCase):
    def test_parses_ps_and_keeps_only_this_run(self):
        a = "/runs/a/Results_Oct20/WorkingDirectory"
        ps = (f"   12 /bin/sh -c iqtree3 -s {a}/Alignments_ids/OG0000001.fa --prefix {a}/Alignments_ids/OG0000001\n"
              f" 7200 iqtree3 -s {a}/Alignments_ids/OG0000000.fa --prefix {a}/Alignments_ids/OG0000000 -quiet\n"
              f"   30 iqtree3 -s {a}/Alignments_ids/OG0000005.fa --prefix {a}/Alignments_ids/OG0000005 -quiet\n"
              "   40 iqtree3 -s /runs/b/Results_Oct20/WorkingDirectory/Alignments_ids/OG0000002.fa -quiet\n"
              "    5 python something\n")

        def fake(cmd, **kw):
            return subprocess.CompletedProcess(cmd, 0, stdout=ps, stderr="")

        self.assertEqual(status.running_iqtree(a, fake), [("OG0000000", 7200), ("OG0000005", 30)])


class TestCountFasta(unittest.TestCase):
    def test_counts_records(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "a.fa"
            fasta(p, 7)
            self.assertEqual(status._count_fasta(p), 7)
            p.write_text("")
            self.assertEqual(status._count_fasta(p), 0)


if __name__ == "__main__":
    unittest.main()


class TestPlanModel(unittest.TestCase):
    def test_giants_use_the_anchor_and_test_halves_them(self):
        from cichlid_orthofinder import status as S
        self.assertAlmostEqual(S.single_seconds(500), S.C * 500 ** S.B)
        self.assertAlmostEqual(S.single_seconds(576 * 4), 12.5 * 3600 * 8)          # (4)^1.5 = 8
        self.assertAlmostEqual(S.single_seconds(5613, "TEST"), S.single_seconds(5613) / 2)
        self.assertLess(S.single_seconds(5613) / 3600 / 24, 20)                      # one core; not the raw curve's months

    def test_simulate_waits_for_threads(self):
        from cichlid_orthofinder import status as S
        # 2 cores: a running 1-core tree with 10 s left; then a 1-core 5 s tree, then a 2-core 4 s tree
        f = S.simulate([("a", 3, 5.0, 1), ("b", 9, 4.0, 2)], [("r", 4, 10.0, 1)], 2)
        self.assertEqual(f["r"], (4, 10.0))
        self.assertEqual(f["a"], (3, 5.0))
        self.assertEqual(f["b"], (9, 14.0))       # waits for both cores (at 10 s), then 4 s
