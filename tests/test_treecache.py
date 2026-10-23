import json
import os
import sys
import tempfile
import threading
import unittest
from pathlib import Path

from cichlid_orthofinder import treecache

FAKE = str(Path(__file__).with_name("fake_iqtree.py"))
ALN = ">0_0\nMKV-L\n>0_1\nMKVAL\n>1_0\nMRVAL\n>2_0\nM-VAL\n"


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.cache, self.work = d / "cache", d / "work"
        self.work.mkdir()
        self.aln = self.work / "OG0000001.fa"
        self.aln.write_text(ALN)
        self.calls = d / "calls.txt"
        self.env = {k: os.environ.get(k) for k in ("FAKE_IQTREE_CALLS", "FAKE_IQTREE_CRASH")}
        os.environ["FAKE_IQTREE_CALLS"] = str(self.calls)
        os.environ.pop("FAKE_IQTREE_CRASH", None)
        os.chmod(FAKE, 0o755)

    def tearDown(self):
        for k, v in self.env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        self.tmp.cleanup()

    def settings(self, seed=1, threads="1", model="MFP", live=None):
        return treecache.Settings(self.cache, seed, threads, FAKE, model=model, live=live)

    def build(self, crash=None, out="OG0000001", **kw):
        if crash:
            os.environ["FAKE_IQTREE_CRASH"] = crash
        else:
            os.environ.pop("FAKE_IQTREE_CRASH", None)
        return treecache.build(self.aln, self.work / out, self.settings(**kw))

    def runs(self):
        return self.calls.read_text().split("\n")[:-1] if self.calls.exists() else []

    def folder(self):
        return next(p.parent for p in self.cache.glob("*/*/aln.fa"))

    def events(self):
        return [json.loads(l) for l in (self.cache / "events.jsonl").read_text().splitlines()]


class TestStates(Case):
    def test_new_then_reused(self):
        ev = self.build()
        self.assertEqual(ev["how"], "built")
        tree = (self.work / "OG0000001.treefile").read_bytes()
        self.assertEqual(len(self.runs()), 1)
        (self.work / "OG0000001.treefile").unlink()
        ev = self.build()
        self.assertEqual(ev["how"], "reused")
        self.assertEqual(len(self.runs()), 1, "a finished tree must not run IQ-TREE again")
        self.assertEqual((self.work / "OG0000001.treefile").read_bytes(), tree)

    def test_killed_during_run_resumes(self):
        with self.assertRaises(SystemExit):
            self.build(crash="checkpoint")
        self.assertFalse((self.work / "OG0000001.treefile").exists(), "no output for an unfinished tree")
        self.assertEqual(treecache.state(self.folder()), "partial")
        ev = self.build()
        self.assertEqual(ev["how"], "resumed")
        self.assertEqual([r.split()[1] for r in self.runs()], ["fresh", "resume"])
        self.assertIn("CHECKPOINT: Resuming", (self.folder() / "run.log").read_text())

    def test_killed_after_iqtree_finished_is_recorded_not_rerun(self):
        self.build()
        (self.folder() / "done.json").unlink()            # as if the wrapper died before recording it
        (self.folder() / "tree.nwk").unlink()
        self.assertEqual(treecache.state(self.folder()), "finished")
        ev = self.build()
        self.assertEqual(ev["how"], "recorded")
        self.assertEqual(len(self.runs()), 1)

    def test_iqtree_killed_at_exit_with_outputs_complete_counts(self):
        ev = self.build(crash="after_tree")
        self.assertEqual(ev["how"], "built")
        self.assertTrue((self.work / "OG0000001.treefile").exists())

    def test_unusable_checkpoint_is_set_aside_and_rebuilt(self):
        with self.assertRaises(SystemExit):
            self.build(crash="checkpoint")
        ev = self.build(crash="bad_resume")
        self.assertEqual(ev["how"], "resumed")             # entered as a resume, finished afresh
        self.assertEqual(len(list(self.folder().glob("failed_*"))), 1)
        self.assertIn("restart_fresh", [e["event"] for e in self.events()])

    def test_failure_raises_and_leaves_no_output(self):
        with self.assertRaises(SystemExit):
            self.build(crash="fail")
        self.assertFalse((self.work / "OG0000001.treefile").exists())
        self.assertEqual(self.events()[-1]["event"], "failed")


class TestKeys(Case):
    def test_seed_change_is_a_miss(self):
        self.build(seed=1)
        self.build(seed=2)
        self.assertEqual(len(self.runs()), 2)
        self.assertEqual(len(list(self.cache.glob("*/*/done.json"))), 2)

    def test_alignment_change_is_a_miss(self):
        self.build()
        self.aln.write_text(ALN.replace("MKVAL", "MKVAI"))
        self.build()
        self.assertEqual(len(self.runs()), 2)

    def test_same_alignment_other_orthogroup_reuses(self):
        self.build(out="OG0000001")
        ev = self.build(out="OG0000999")
        self.assertEqual(ev["how"], "reused")
        self.assertEqual((self.work / "OG0000001.treefile").read_bytes(),
                         (self.work / "OG0000999.treefile").read_bytes())

    def test_threads_by_size_are_in_the_key(self):
        self.assertEqual(treecache.parse_threads("1;500:8;2000:16"), [(0, 1), (500, 8), (2000, 16)])
        steps = treecache.parse_threads("1;500:8;2000:16")
        self.assertEqual([treecache.threads_for(n, steps) for n in (3, 499, 500, 1999, 5613)], [1, 1, 8, 8, 16])
        steps = treecache.parse_threads("1;600:2;1500:8;3000:16")     # the run's default
        self.assertEqual([treecache.threads_for(n, steps) for n in (599, 600, 1499, 1500, 2999, 3000, 5613)],
                         [1, 2, 2, 8, 8, 16, 16])
        self.build(threads="1")
        self.build(threads="1;3:4")                          # 4 sequences -> 4 threads: a different key
        self.assertEqual(len(self.runs()), 2)


class TestConcurrency(Case):
    def test_two_calls_at_once_run_iqtree_once(self):
        results, errors = [], []

        def go(name):
            try:
                results.append(treecache.build(self.aln, self.work / name, self.settings()))
            except BaseException as e:
                errors.append(e)

        ts = [threading.Thread(target=go, args=(n,)) for n in ("OG0000001", "OG0000002")]
        [t.start() for t in ts]
        [t.join() for t in ts]
        self.assertEqual(errors, [])
        self.assertEqual(sorted(r["how"] for r in results), ["built", "reused"])
        self.assertEqual(len(self.runs()), 1)


class TestHelpers(unittest.TestCase):
    def test_alignment_shape(self):
        self.assertEqual(treecache.alignment_shape(ALN.encode()), (4, 5))
        self.assertEqual(treecache.alignment_shape(b">a\nMK\nVL\n>b\nMKVL\n"), (2, 4))

    def test_missing_cache_setting_fails_clearly(self):
        with self.assertRaises(SystemExit) as e:
            treecache.Settings.from_env({})
        self.assertIn("CICHLID_TREE_CACHE", str(e.exception))

    def test_iqtree_args_have_model_and_seed(self):
        a = treecache.iqtree_args("x.fa", "p", 7, 8)
        self.assertEqual(a[a.index("-m") + 1], "MFP")
        self.assertEqual(a[a.index("-seed") + 1], "7")
        self.assertEqual(a[a.index("-T") + 1], "8")


class TestBudget(unittest.TestCase):
    def test_waits_while_pool_is_full_and_reclaims_dead_reservations(self):
        with tempfile.TemporaryDirectory() as d:
            cache = Path(d)
            order = []
            with treecache.reserved(cache, 3, 4, poll=0.05):
                def second():
                    with treecache.reserved(cache, 2, 4, poll=0.05):
                        order.append("second got 2")
                t = threading.Thread(target=second)
                t.start()
                import time
                time.sleep(0.3)
                order.append("first releases")
            t.join(5)
            self.assertEqual(order, ["first releases", "second got 2"])
            # a reservation left by a process that no longer exists is reclaimed
            (cache / "budget" / "999999999-1.res").write_text("4")
            with treecache.reserved(cache, 4, 4, poll=0.05):
                pass
            self.assertFalse((cache / "budget" / "999999999-1.res").exists())

    def test_request_larger_than_pool_is_capped(self):
        with tempfile.TemporaryDirectory() as d:
            with treecache.reserved(Path(d), 16, 4, poll=0.05) as got:
                self.assertEqual(got, 4)



class TestModelRule(Case):
    """The model setting by size, and the live settings file that changes it while a run goes on."""

    def test_rule_parsing(self):
        steps = treecache.parse_models("MFP;1000:TEST")
        self.assertEqual(steps, [(0, "MFP"), (1000, "TEST")])
        self.assertEqual([treecache.model_for(n, steps) for n in (3, 999, 1000, 5613)], ["MFP", "MFP", "TEST", "TEST"])
        self.assertEqual(treecache.parse_models("MFP"), [(0, "MFP")])
        self.assertEqual(treecache.parse_models("LG+G4"), [(0, "LG+G4")])
        for bad in ("MFP;TEST", "MFP;1000:", "MFP;1000:TE ST", "MFP;x:TEST", "MFP;100:TEST;100:MFP"):
            with self.assertRaises(ValueError, msg=bad):
                treecache.parse_models(bad)

    def test_mfp_keeps_the_key_it_had(self):
        # trees built with MFP keep their cache key, so a cache from before the model rule is still used
        data = self.aln.read_bytes()
        old = {"program": "IQ-TREE", "version": "9.9.9", "model": "MFP", "seed": 1, "threads": 1}
        self.build(model="MFP;1000:TEST")
        self.assertEqual(self.folder().name, treecache.cache_key(old, data))

    def test_model_by_size_reaches_iqtree(self):
        self.build(model="MFP;4:TEST")                   # the test alignment has 4 sequences
        self.assertEqual(self.runs()[-1].split()[2], "TEST")
        ev = self.events()[-1]
        self.assertEqual(ev["model"], "TEST")
        done = json.loads((self.folder() / "done.json").read_text())
        self.assertEqual(done["params"]["model"], "TEST")

    def test_live_file_changes_trees_that_start_later(self):
        live = Path(self.tmp.name) / "tree_settings.json"
        live.write_text(json.dumps({"model": "MFP"}))
        self.build(model="MFP;4:TEST", live=live)       # the file wins over the launch setting
        self.assertEqual(self.runs()[-1].split()[2], "MFP")
        live.write_text(json.dumps({"model": "MFP;4:TEST"}))     # changed while the run goes on
        self.aln.write_text(ALN.replace("MKVAL", "MKVAV"))        # the next tree: another alignment
        self.build(out="OG0000002", model="MFP", live=live)
        self.assertEqual(self.runs()[-1].split()[2], "TEST")
        self.assertEqual(len(self.runs()), 2)

    def test_finished_tree_is_not_rebuilt_under_the_same_model(self):
        live = Path(self.tmp.name) / "tree_settings.json"
        live.write_text(json.dumps({"model": "MFP;4:TEST"}))
        self.build(live=live)
        self.build(live=live)
        self.assertEqual(len(self.runs()), 1)

    def test_unreadable_live_file_falls_back_and_is_recorded(self):
        live = Path(self.tmp.name) / "tree_settings.json"
        live.write_text("{not json")
        self.build(model="MFP;4:TEST", live=live)
        self.assertEqual(self.runs()[-1].split()[2], "TEST")
        self.assertIn("settings_unreadable", [e["event"] for e in self.events()])
        live.write_text(json.dumps({"model": "MFP;1000:TE ST"}))   # a bad rule: also falls back
        self.aln.write_text(ALN.replace("MKVAL", "MKVAV"))
        self.build(out="OG0000002", model="MFP;4:TEST", live=live)
        self.assertEqual(self.runs()[-1].split()[2], "TEST")


if __name__ == "__main__":
    unittest.main()
