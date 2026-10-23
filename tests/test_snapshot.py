import tempfile
import unittest
import unittest.mock
from pathlib import Path

from cichlid_orthofinder import snapshot

LOG = """2026-10-22 : Started OrthoFinder version 3.1.5
Command Line: orthofinder -f /orig/input -o /orig/out
WorkingDirectory_Base: /orig/out/Results_Oct22/WorkingDirectory/
FN_Orthogroups: /orig/out/Results_Oct22/WorkingDirectory/clusters_OrthoFinder_I1.2.txt_id_pairs.txt
WorkingDirectory_Trees: /orig/out/Results_Oct22/WorkingDirectory/OrthoFinder/Results_Oct22/WorkingDirectory/
"""


class TestLog(unittest.TestCase):
    def test_rewrite_points_all_three_paths_at_the_copy(self):
        with tempfile.TemporaryDirectory() as d:
            src, dst = Path(d) / "src", Path(d) / "dst"
            (dst / "WD").mkdir(parents=True), src.mkdir()
            (dst / "WD" / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text("")
            (src / "Log.txt").write_text(LOG)
            self.assertTrue(snapshot.log_paths(src)["WorkingDirectory_Trees"].endswith("Results_Oct22/WorkingDirectory/"))
            snapshot._rewrite_log(src, dst, dst / "WD")
            q = snapshot.log_paths(dst)
            self.assertEqual(q["WorkingDirectory_Base"], f"{dst / 'WD'}/")
            self.assertEqual(q["WorkingDirectory_Trees"], f"{dst / 'WD'}/")
            self.assertEqual(q["FN_Orthogroups"], str(dst / "WD" / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt"))
            self.assertNotIn("/orig/out/Results_Oct22/WorkingDirectory/", (dst / "Log.txt").read_text().split("Command Line")[1].split("\n", 1)[1])

    def test_log_of_a_running_run_has_no_trees_path_and_it_is_added(self):
        with tempfile.TemporaryDirectory() as d:
            src, dst = Path(d) / "src", Path(d) / "dst"
            (dst / "WD").mkdir(parents=True), src.mkdir()
            (dst / "WD" / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text("")
            (src / "Log.txt").write_text("\n".join(l for l in LOG.splitlines() if not l.startswith("WorkingDirectory_Trees")))
            self.assertNotIn("WorkingDirectory_Trees", snapshot.log_paths(src))
            snapshot._rewrite_log(src, dst, dst / "WD")
            self.assertEqual(snapshot.log_paths(dst)["WorkingDirectory_Trees"], f"{dst / 'WD'}/")

    def test_no_clusters_file_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "Log.txt").write_text(LOG)
            with self.assertRaises(RuntimeError):
                snapshot._rewrite_log(d, d, Path(d))


def fasta(path, n):
    path.write_text("".join(f">0_{i}\nMKV\n" for i in range(n)))


class TestTakeWhileRunning(unittest.TestCase):
    """A running (and a resumed) run: Log.txt has no WorkingDirectory_Trees yet; trees are read from the
    current Results folder; the species files from WorkingDirectory_Base."""

    def make(self, d, resumed):
        d = Path(d)
        first = d / "out" / "Results_Oct22"
        base = first / "WorkingDirectory"
        base.mkdir(parents=True)
        for f in ("SpeciesIDs.txt", "SequenceIDs.txt", "Species0.fa", "Blast0_0.txt.gz"):
            (base / f).write_text(f)
        current_results = (base / "OrthoFinder" / "Results_Oct22") if resumed else first
        cur = current_results / "WorkingDirectory"
        for sub in ("Sequences_ids", "Trees_ids", "Alignments_ids"):
            (cur / sub).mkdir(parents=True, exist_ok=True)
        (cur / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text("")
        for og, n, done in (("OG0000000", 50, False), ("OG0000001", 8, True), ("OG0000002", 5, True)):
            fasta(cur / "Sequences_ids" / f"{og}.fa", n)
            (cur / "Alignments_ids" / f"{og}.fa").write_text(">a\nM\n")
            if done:
                (cur / "Trees_ids" / f"{og}.txt").write_text("(a,b,c);")
        (current_results / "Log.txt").write_text(f"Command Line: x\nWorkingDirectory_Base: {base}/\n")
        return current_results

    def check(self, resumed):
        with tempfile.TemporaryDirectory() as d:
            results = self.make(d, resumed)

            def fake_orthofinder(cmd, **kw):            # stands in for `orthofinder -fgt <stage>`
                stage = Path(cmd[cmd.index("-fgt") + 1])
                q = snapshot.log_paths(stage)
                assert set(q) == set(snapshot.LOG_KEYS), q
                (stage.parent / "Results_Oct23" / "Orthogroups").mkdir(parents=True)
                return unittest.mock.Mock(returncode=0)

            made, manifest = snapshot.take(results, Path(d) / "snap", runner=fake_orthofinder, say=lambda m: None)
            wd = Path(d) / "snap" / "Results_snapshot" / "WorkingDirectory"
            self.assertEqual(manifest["trees_included"], 2)
            self.assertEqual([e["cluster"] for e in manifest["excluded"]], ["OG0000000"])
            self.assertTrue((wd / "SpeciesIDs.txt").is_file())               # from the first working folder
            self.assertFalse((wd / "Blast0_0.txt.gz").exists())              # the search is not copied
            self.assertEqual(sorted(p.name for p in (wd / "Trees_ids").iterdir()), ["OG0000001.txt", "OG0000002.txt"])
            self.assertEqual(made.name, "Results_Oct23")

    def test_running_fresh_run(self):
        self.check(resumed=False)

    def test_running_resumed_run(self):
        self.check(resumed=True)


class TestMapping(unittest.TestCase):
    def results(self, root, ogs, trees):
        (root / "Orthogroups").mkdir(parents=True)
        (root / "Resolved_Gene_Trees").mkdir()
        (root / "Orthogroups" / "Orthogroups.txt").write_text("".join(f"{k}: {' '.join(v)}\n" for k, v in ogs.items()))
        (root / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text(
            "".join(f"{k}: {t}\n" for k, t in trees.items()))
        return root

    def test_maps_by_genes_and_compares_trees(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            snap = self.results(d / "s", {"OG0000000": ["A_1", "B_1"], "OG0000001": ["A_2", "B_2", "C_2"],
                                          "OG0000002": ["A_9"]},
                                {"OG0000000": "(A_A_1,B_B_1);", "OG0000001": "(A_A_2,(B_B_2,C_C_2));"})
            final = self.results(d / "f", {"OG0000000": ["X_1", "Y_1", "Z_1", "W_1"], "OG0000001": ["A_2", "B_2", "C_2"],
                                           "OG0000002": ["A_1", "B_1"]},
                                 {"OG0000001": "(A_A_2,(B_B_2,C_C_2));", "OG0000002": "(A_A_1,B_B_1);"})
            s = snapshot.mapping(snap, final, d / "m.tsv")
            self.assertEqual(s, {"orthogroups": 2, "matched": 2, "same_tree": 2})   # the one-gene group is skipped
            rows = (d / "m.tsv").read_text().splitlines()[1:]
            self.assertEqual(rows, ["OG0000000\tOG0000002\t2\tyes", "OG0000001\tOG0000001\t3\tyes"])


if __name__ == "__main__":
    unittest.main()
