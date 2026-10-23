import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from cichlid_orthofinder import run, snapshot


def fasta(path, n):
    path.write_text("".join(f">0_{i}\nMKV\n" for i in range(n)))


def results_with(root, sizes, done):
    r = Path(root) / "output" / "Results_Oct22"
    for d in ("Sequences_ids", "Trees_ids"):
        (r / "WorkingDirectory" / d).mkdir(parents=True)
    for og, n in sizes.items():
        fasta(r / "WorkingDirectory" / "Sequences_ids" / f"{og}.fa", n)
    for og in done:
        (r / "WorkingDirectory" / "Trees_ids" / f"{og}.txt").write_text("(a,b,c);")
    return r


SIZES = {"OG0000000": 1500, "OG0000001": 800, "OG0000002": 400, "OG0000003": 20, "OG0000004": 5, "OG0000005": 2}


def finish_tree(r, og):
    (r / "WorkingDirectory/Trees_ids" / f"{og}.txt").write_text("(a,b,c);")


class TestReady(unittest.TestCase):
    def test_up_to_and_including_the_size(self):
        with tempfile.TemporaryDirectory() as d:
            r = results_with(d, SIZES, ["OG0000003", "OG0000004"])
            self.assertFalse(snapshot.ready(r, 500))            # the 400 still running
            finish_tree(r, "OG0000002")
            self.assertTrue(snapshot.ready(r, 500))             # all up to 500 done, larger ones running
            self.assertTrue(snapshot.ready(r, 400))             # inclusive: the 400 counts as "up to 400"
            self.assertFalse(snapshot.ready(r, 1000))           # the 800 still running
            finish_tree(r, "OG0000001")
            finish_tree(r, "OG0000000")
            self.assertFalse(snapshot.ready(r, 1000))           # nothing larger left: the final results come next

    def test_not_before_the_tree_stage(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertFalse(snapshot.ready(Path(d), 500))


def fake_take(results, dest, *a, **k):
    return Path(dest) / "Results_X", {"taken": str(dest), "clusters_excluded": 1, "genes_excluded": 9}


class TestWatcher(unittest.TestCase):
    def watcher(self, d, state=None):
        return snapshot.Watcher(Path(d), Path(d) / "output", [500, 1000], None, state=state, say=lambda m: None)

    def test_one_snapshot_per_size_in_order_then_maps(self):
        with tempfile.TemporaryDirectory() as d:
            r = results_with(d, SIZES, ["OG0000002", "OG0000003", "OG0000004"])
            state = {}
            w = self.watcher(d, state)
            with unittest.mock.patch.object(snapshot, "take", side_effect=fake_take) as take:
                self.assertTrue(w.poll())                        # up to 500
                self.assertFalse(w.poll())                       # not again; 1000 not ready
                finish_tree(r, "OG0000001")
                self.assertTrue(w.poll())                        # up to 1000
                finish_tree(r, "OG0000000")
                self.assertFalse(w.poll())                       # done
                self.assertEqual(take.call_count, 2)
            self.assertEqual([s["size"] for s in state["snapshots"]], [500, 1000])
            with unittest.mock.patch.object(snapshot, "mapping",
                                            return_value={"orthogroups": 3, "matched": 3, "same_tree": 3}) as m:
                w.finish("/final")
                self.assertEqual(m.call_count, 2)
            self.assertTrue(all(s["mapping_same_tree"] == 3 for s in state["snapshots"]))

    def test_resume_does_not_repeat_and_skips_an_overtaken_size(self):
        with tempfile.TemporaryDirectory() as d:
            # a resume finds everything up to 1000 done: only the 1000 snapshot is taken
            results_with(d, SIZES, ["OG0000001", "OG0000002", "OG0000003", "OG0000004"])
            state = {}
            w = self.watcher(d, state)
            with unittest.mock.patch.object(snapshot, "take", side_effect=fake_take) as take:
                self.assertTrue(w.poll())
                self.assertEqual(take.call_count, 1)
            self.assertEqual(state["snapshots"][0]["size"], 1000)
            self.assertEqual(state["snapshots"][0]["skipped_smaller"], [500])
            # a new watcher with the saved state (the next attempt) takes nothing
            w2 = self.watcher(d, state)
            with unittest.mock.patch.object(snapshot, "take") as take:
                self.assertFalse(w2.poll())
                take.assert_not_called()


class TestRunWatched(unittest.TestCase):
    def test_polls_while_running_and_survives_a_failing_snapshot(self):
        with tempfile.TemporaryDirectory() as d:
            calls, said = [], []

            class W:
                def poll(self):
                    calls.append(1)
                    raise RuntimeError("disk full")

            cmd = [sys.executable, "-c", "import time; time.sleep(0.6); print('Done writing files')"]
            rc = run.run_watched(cmd, None, Path(d) / "log.txt", W(), poll_seconds=0.1, say=said.append)
            self.assertEqual(rc, 0)
            self.assertGreaterEqual(len(calls), 2)
            self.assertTrue(any("run continues" in s for s in said))
            self.assertTrue(run.finished_ok(rc, Path(d) / "log.txt"))

    def test_exit_zero_with_a_traceback_is_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "log.txt"
            log.write_text("Traceback (most recent call last):\n  ...\nOrthoFinder finished in 4.4s\n")
            self.assertFalse(run.finished_ok(0, log))


class TestMissingTrees(unittest.TestCase):
    def test_finished_but_missing_trees_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            r = results_with(d, SIZES, ["OG0000002", "OG0000003", "OG0000004"])
            self.assertEqual(run.missing_trees(r), [("OG0000000", 1500), ("OG0000001", 800)])
            finish_tree(r, "OG0000000"); finish_tree(r, "OG0000001")
            self.assertEqual(run.missing_trees(r), [])             # OG0000005 has 2 sequences: no tree expected


class TestCompleteness(unittest.TestCase):
    def test_success_needs_every_tree_and_no_failed_program(self):
        with tempfile.TemporaryDirectory() as d:
            r = results_with(d, {"OG0000000": 10, "OG0000001": 5, "OG0000002": 2}, ["OG0000001"])
            log = Path(d) / "log.txt"
            log.write_text("Done writing files\nOrthoFinder finished in 23s\n")
            self.assertEqual(run.missing_trees(r), [("OG0000000", 10)])   # the 2-sequence cluster needs none
            self.assertFalse(run.finished_ok(0, log, r))                    # a tree is missing
            (r / "WorkingDirectory/Trees_ids/OG0000000.txt").write_text("(a,b,c);")
            self.assertTrue(run.finished_ok(0, log, r))
            log.write_text("ERROR: external program returned code -9\nDone writing files\n")
            self.assertFalse(run.finished_ok(0, log, r))                    # a failed external program


if __name__ == "__main__":
    unittest.main()
