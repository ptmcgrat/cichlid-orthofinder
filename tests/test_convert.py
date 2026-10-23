"""The converter, on a fixture cut from pilot 2's real OrthoFinder output (tests/fixtures/pilot2/: five whole
preliminary clusters, the ten orthogroups cut from them and four unassigned genes, 55 proteins)."""

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from cichlid_data import schemas

from cichlid_orthofinder import convert as C

FIXTURE = Path(__file__).with_name("fixtures") / "pilot2"
SPECIES = ["Mzeb", "Acal", "Gacu", "Drer", "Locu"]          # any order: the catalog's order in the real run


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        shutil.copytree(FIXTURE, self.root / "in")
        self.results = self.root / "in" / "Results_Oct20"
        self.working = self.results / "WorkingDirectory"

    def tearDown(self):
        self.tmp.cleanup()

    def out(self, name="a"):
        return self.root / name / "Results" / "v1" / "OrthoFinder"

    def convert(self, name="a", **kw):
        proteomes = {p: self.root / "in" / "proteomes" / f"{p}.faa" for p in SPECIES}
        return C.convert(self.results, self.out(name), SPECIES, proteomes, "3.1.5", say=lambda *a: None, **kw)


class TestOutputs(Case):
    def test_every_gene_in_its_assign_file_in_proteome_order(self):
        self.convert()
        for p in SPECIES:
            rows = schemas.ASSIGN.read(self.out() / "assign" / f"{p}.assign.tsv")
            proteome = [n for n, _ in C.read_fasta(self.root / "in" / "proteomes" / f"{p}.faa")]
            self.assertEqual([r["gene_id"] for r in rows], proteome)
            self.assertEqual(schemas.ASSIGN.validate_file(self.out() / "assign" / f"{p}.assign.tsv"), [])
            for r in rows:
                self.assertTrue(r["protein_length"])
                self.assertEqual({k for k, v in r.items() if v}, {"gene_id", "protein_length"} | (
                    {"best_family"} if r["best_family"] else set()))

    def test_gene_in_no_orthogroup_has_an_empty_family(self):
        self.convert()
        rows = {r["gene_id"]: r for r in schemas.ASSIGN.read(self.out() / "assign" / "Gacu.assign.tsv")}
        self.assertEqual(rows["Gacu_120815963"]["best_family"], "")          # unassigned in OrthoFinder
        acal = {r["gene_id"]: r for r in schemas.ASSIGN.read(self.out() / "assign" / "Acal.assign.tsv")}
        self.assertEqual(acal["Acal_113030444"]["best_family"], "OG0000339")  # a 3-gene orthogroup: no folder

    def test_families_validate_with_our_ids_only(self):
        summary = self.convert()
        self.assertEqual(sorted(p.name for p in (self.out() / "families").iterdir()),
                         ["OG0000005", "OG0000012", "OG0000065", "OG0000071", "OG0000107", "OG0000112"])
        self.assertEqual(schemas.validate_families(self.out() / "families", prefixes=set(SPECIES)), [])
        self.assertEqual(summary["counts"]["orthogroups_below_min_size"], 4)
        for fam in (self.out() / "families").iterdir():
            leaves = [l.name for l in C.parse_newick((fam / "tree.nwk").read_text()).leaves()]
            self.assertFalse([l for l in leaves if C.INTERNAL_RE.match(l)])
            self.assertTrue(all(l.split("_", 1)[0] in SPECIES and l.count("_") == 1 for l in leaves))

    def test_done_json_is_written_last_and_has_only_the_contract_keys(self):
        self.convert()
        for fam in (self.out() / "families").iterdir():
            done = json.loads((fam / "done.json").read_text())
            self.assertEqual(set(done), {"release", "rooted", "tree_mode", "inputs_hash"})
            self.assertEqual((done["release"], done["rooted"], done["tree_mode"]), ("v1", False, "denovo"))
            self.assertGreaterEqual((fam / "done.json").stat().st_mtime_ns, (fam / "tree.nwk").stat().st_mtime_ns)
        self.assertFalse(list((self.out() / "families").glob(".*")), "no half-written family is left")

    def test_unrooted_raw_and_rooted_resolved(self):
        self.convert("raw")
        self.convert("res", tree_set="resolved")
        for fam in (self.out("raw") / "families").iterdir():
            self.assertEqual(len(C.parse_newick((fam / "tree.nwk").read_text()).children), 3)
        for fam in (self.out("res") / "families").iterdir():
            done = json.loads((fam / "done.json").read_text())
            self.assertEqual(len(C.parse_newick((fam / "tree.nwk").read_text()).children), 2)
            self.assertEqual((done["rooted"], done["rooting"]), (True, "orthofinder species tree"))
        self.assertEqual(schemas.validate_families(self.out("res") / "families", prefixes=set(SPECIES)), [])

    def test_raw_cluster_alignment_matches_the_tree(self):
        self.convert(alignment_source="raw_cluster")
        self.assertEqual(schemas.validate_families(self.out() / "families", prefixes=set(SPECIES)), [])
        rows = C.read_fasta(self.out() / "families" / "OG0000107" / "combined.aln")
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(any(s[i] != "-" for _, s in rows) for i in range(len(rows[0][1]))))

    def test_converting_twice_gives_identical_bytes(self):
        for name in ("a", "b"):
            self.convert(name)
        files_a = sorted(p.relative_to(self.out("a")) for p in self.out("a").rglob("*") if p.is_file())
        files_b = sorted(p.relative_to(self.out("b")) for p in self.out("b").rglob("*") if p.is_file())
        self.assertEqual(files_a, files_b)
        for rel in files_a:
            self.assertEqual((self.out("a") / rel).read_bytes(), (self.out("b") / rel).read_bytes(), str(rel))

    def test_summary_tables(self):
        self.convert()
        split = (self.out() / "native" / "summary" / "split_clusters.tsv").read_text().splitlines()
        self.assertEqual(split[0].split("\t")[:2], ["cluster", "orthogroup"])
        self.assertEqual(len(split) - 1, 8)                    # the pieces of clusters 0, 5 and 52
        models = (self.out() / "native" / "summary" / "cluster_models.tsv").read_text().splitlines()
        self.assertEqual(models[1:], ["OG0000000\t20\tMFP", "OG0000005\t15\tMFP", "OG0000052\t8\tMFP",
                                      "OG0000079\t6\tMFP", "OG0000316\t3\tMFP"])


class TestRefusals(Case):
    def test_leaf_that_cannot_be_mapped(self):
        f = self.working / "Trees_ids" / "OG0000079.txt"
        f.write_text(f.read_text().replace("(", "(9_999:0.1,", 1))
        with self.assertRaisesRegex(C.ConversionError, "9_999.*not in SequenceIDs"):
            self.convert()

    def test_label_without_the_doubled_prefix(self):
        self.assertEqual(C.strip_species("Acal_Acal_113012997", SPECIES), "Acal_113012997")
        with self.assertRaises(C.ConversionError):
            C.strip_species("Acal_113012997", SPECIES)
        f = self.results / "MultipleSequenceAlignments" / "OG0000071.fa"
        f.write_text(f.read_text().replace("Acal_Acal_", "Acal_", 1))
        with self.assertRaisesRegex(C.ConversionError, "OG0000071.fa: row 'Acal_113031360' is not OrthoFinder's"):
            self.convert()

    def test_one_row_none_alignment(self):
        self.assertEqual([n for n, _ in C.read_fasta(self.results / "MultipleSequenceAlignments" / "OG0000339.fa")],
                         ["None"])
        genes = ["Acal_113030444", "Drer_554064", "Gacu_120817097"]
        g2i = {g: i for i, g in C.read_sequence_ids(self.working).items()}
        with self.assertRaisesRegex(C.ConversionError, "OG0000339.fa: row 'None'"):
            C.alignment(self.results, "orthogroup", "OG0000339", 316, genes, SPECIES, g2i)
        rows = C.alignment(self.results, "raw_cluster", "OG0000339", 316, genes, SPECIES, g2i)
        self.assertEqual([n for n, _ in rows], sorted(genes))      # the cluster's own alignment has the rows

    def test_small_orthogroup_without_an_alignment_file(self):
        # in pilot 2, 2- and 3-gene pieces of split clusters have no MultipleSequenceAlignments file
        self.assertFalse((self.results / "MultipleSequenceAlignments" / "OG0000337.fa").exists())
        with self.assertRaisesRegex(C.ConversionError, "OG0000337: no MultipleSequenceAlignments"):
            self.convert(min_size=3)

    def test_duplicate_rows_and_unequal_lengths(self):
        with self.assertRaisesRegex(C.ConversionError, "duplicate rows"):
            C.check_alignment([("a", "MK"), ("a", "MK")], ["a"], "x")
        with self.assertRaisesRegex(C.ConversionError, "different lengths"):
            C.check_alignment([("a", "MK"), ("b", "M")], ["a", "b"], "x")

    def test_orthogroup_with_no_tree_gets_no_folder(self):
        (self.working / "Trees_ids" / "OG0000079.txt").unlink()
        summary = self.convert()
        self.assertFalse((self.out() / "families" / "OG0000071").exists())
        self.assertEqual(summary["counts"]["orthogroups_without_tree"], 1)
        rows = {r["gene_id"]: r for r in schemas.ASSIGN.read(self.out() / "assign" / "Acal.assign.tsv")}
        self.assertEqual(rows["Acal_113031360"]["best_family"], "OG0000071")   # still assigned

    def test_species_or_proteins_that_do_not_match(self):
        with self.assertRaisesRegex(C.ConversionError, "SpeciesIDs.txt has"):
            C.convert(self.results, self.out(), SPECIES[:4], {}, "3.1.5", say=lambda *a: None)
        f = self.root / "in" / "proteomes" / "Locu.faa"
        f.write_text(f.read_text() + ">Locu_1\nMK\n")
        with self.assertRaisesRegex(C.ConversionError, "sequences, the proteomes"):
            self.convert()

    def test_orthogroup_not_in_exactly_one_cluster(self):
        f = self.working / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt"
        text = f.read_text()
        first = text.split("begin", 1)[1].split()[1]           # move one gene of cluster 0 to cluster 5
        text = text.replace(f" {first} ", " ", 1).replace("\n5 ", f"\n5 {first} ", 1)
        f.write_text(text)
        with self.assertRaisesRegex(C.ConversionError, "exactly one preliminary cluster"):
            self.convert()

    def test_refuses_family_folders_of_another_conversion(self):
        self.convert()
        stale = self.out() / "families" / "OG0000107"
        with self.assertRaisesRegex(C.ConversionError, "2 family folders are not part of this conversion"):
            self.convert(min_size=6)                          # OG0000107 and OG0000112 (5 genes) would not be now
        self.assertTrue((stale / "done.json").is_file(), "nothing is deleted")


class TestRerun(Case):
    def files(self, name="a"):
        return {p.relative_to(self.out(name)): p.read_bytes() for p in self.out(name).rglob("*") if p.is_file()}

    def test_rerun_reuses_finished_families_and_gives_the_same_files(self):
        first = self.convert()
        before = self.files()
        mtimes = {f.name: (f / "done.json").stat().st_mtime_ns for f in (self.out() / "families").iterdir()}
        again = self.convert()
        self.assertEqual(first["families_by_how"], {"built": 6})
        self.assertEqual(again["families_by_how"], {"reused": 6})
        self.assertEqual(self.files(), before)
        self.assertEqual({f.name: (f / "done.json").stat().st_mtime_ns for f in (self.out() / "families").iterdir()},
                         mtimes, "a reused family is not rewritten")

    def test_unfinished_family_is_built_again(self):
        self.convert()
        before = self.files()
        (self.out() / "families" / "OG0000071" / "done.json").unlink()       # as if stopped before done.json
        again = self.convert()
        self.assertEqual(again["families_by_how"], {"reused": 5, "rebuilt": 1})
        self.assertEqual(self.files(), before)

    def test_family_with_changed_inputs_is_rebuilt(self):
        self.convert()
        (self.out() / "families" / "OG0000071" / "tree.nwk").write_text("(a,b,c);\n")
        f = self.out() / "families" / "OG0000065" / "done.json"
        done = json.loads(f.read_text())
        done["inputs_hash"] = "0" * 64
        f.write_text(json.dumps(done))
        again = self.convert()
        self.assertEqual(again["families_by_how"], {"reused": 4, "rebuilt": 2})
        self.assertEqual(schemas.validate_families(self.out() / "families", prefixes=set(SPECIES)), [])
        self.assertNotEqual((self.out() / "families" / "OG0000071" / "tree.nwk").read_text(), "(a,b,c);\n")
        self.assertFalse(list((self.out() / "families").glob(".*")))


class TestSummaryAndPaths(Case):
    def test_counts_and_skipped_orthogroups_in_native_summary(self):
        self.convert()
        summary = self.out() / "native" / "summary"
        counts = dict(l.split("\t") for l in (summary / "conversion_counts.tsv").read_text().splitlines()[1:])
        self.assertEqual((counts["families_written"], counts["orthogroups_below_min_size"],
                          counts["orthogroups_without_tree"], counts["min_family_size"]), ("6", "4", "0", "4"))
        skipped = (summary / "skipped_orthogroups.tsv").read_text().splitlines()
        self.assertEqual(skipped[0], "orthogroup\tgenes\treason")
        self.assertEqual([l.split("\t")[0] for l in skipped[1:]], ["OG0000337", "OG0000339", "OG0000345", "OG0000346"])
        self.assertTrue(all(l.endswith("fewer than 4 genes") for l in skipped[1:]))

    def test_no_tree_is_listed_as_skipped(self):
        (self.working / "Trees_ids" / "OG0000079.txt").unlink()
        self.convert()
        skipped = (self.out() / "native" / "summary" / "skipped_orthogroups.tsv").read_text()
        self.assertIn("OG0000071\t6\tno gene tree", skipped)

    def test_paths_come_from_the_results_layout(self):
        class Where:                                          # stands in for a coordinator DataPath
            def __init__(self, rel):
                self.local = Path("/R") / rel

        class FakeResultsLayout:
            def assign(self, p): return Where(f"assign/{p}.assign.tsv")
            def family(self, f): return Where(f"families/{f}")
            def section(self, s): return Where(s)
            def native(self): return Where("native")
            def run_record(self): return Where("run.json")

        paths = C.OutputPaths.of(FakeResultsLayout())
        self.assertEqual(paths.layout.__class__.__name__, "FakeResultsLayout")
        self.assertEqual((paths.assign("Mzeb"), paths.family("OG1"), paths.families(), paths.native(),
                          paths.run_record()),
                         (Path("/R/assign/Mzeb.assign.tsv"), Path("/R/families/OG1"), Path("/R/families"),
                          Path("/R/native"), Path("/R/run.json")))
        scratch = C.OutputPaths.of(self.root / "x")
        self.assertEqual(scratch.family("OG1"), self.root / "x" / "families" / "OG1")
        with self.assertRaises(ValueError):
            C.OutputPaths()


class TestTrees(unittest.TestCase):
    TREE = "((A:1,B:2):0.5,(C:1,(D:3,E:1):0.25):0.75,F:4);"

    def test_newick_round_trip_keeps_length_strings(self):
        t = "(0_1:0.0000010782,(0_2:0.0092958990,1_3:1e-06):0.5,2_4:2.1428e-06);"
        self.assertEqual(C.write_newick(C.parse_newick(t)), t)

    def test_pruning_keeps_pairwise_path_lengths(self):
        before = C.path_lengths(C.parse_newick(self.TREE))
        for keep in ({"A", "C", "E", "F"}, {"B", "D", "E"}, {"A", "B", "C", "D"}):
            root = C.unroot(C.prune(C.parse_newick(self.TREE), keep))
            after = C.path_lengths(root)
            self.assertEqual({l.name for l in root.leaves()}, keep)
            for pair, d in after.items():
                self.assertAlmostEqual(d, before[pair], places=9, msg=f"{keep} {pair}")

    def test_no_degree_two_nodes_after_pruning(self):
        root = C.prune(C.parse_newick(self.TREE), {"A", "C", "D", "F"})
        for node in root.iter():
            if node is not root:
                self.assertNotEqual(len(node.children), 1)

    def test_two_child_root_after_pruning_is_unrooted(self):
        root = C.prune(C.parse_newick(self.TREE), {"C", "D", "E", "A"})
        self.assertEqual(len(root.children), 2)                # pruning left a two-child root
        root = C.unroot(root)
        self.assertEqual(len(root.children), 3)
        self.assertEqual(sorted(l.name for l in root.leaves()), ["A", "C", "D", "E"])

    def test_sums_use_ten_significant_digits(self):
        root = C.prune(C.parse_newick("((A:0.1,B:0.2):0.30000000004,C:0.5,D:0.7);"), {"A", "C", "D"})
        self.assertEqual(C.write_newick(root), "(A:0.4,C:0.5,D:0.7);")


class TestRunRecord(Case):
    def test_run_record_validates(self):
        summary = self.convert()
        settings = self.root / "run_settings.json"
        settings.write_text(json.dumps({"attempts": [
            {"command": ["cichlid-orthofinder-of", "-f", "in"], "started": "2026-10-22T18:24:05-04:00",
             "finished": "2026-10-23T12:42:22-04:00", "versions": {"orthofinder": "3.1.5", "iqtree3": "3.1.3"},
             "settings": {"seed": 1, "model": {"value": "full", "rule": "MFP"}}}]}))
        rec = C.run_record(settings, summary, "0.14.1", "0" * 32, "f" * 32, SPECIES)
        schemas.write_run(self.out() / "run.json", **rec)
        self.assertEqual(schemas.validate_run(self.out() / "run.json"), [])
        self.assertEqual(rec["method_version"], "3.1.5")
        self.assertEqual(rec["parameters"]["seed"], 1)
        self.assertEqual(rec["parameters"]["tree_set"], "raw_pruned")
        self.assertEqual(rec["parameters"]["counts"]["families_written"], 6)


if __name__ == "__main__":
    unittest.main()


class TestRunRecordModelHistory(Case):
    def test_model_history_comes_from_tree_settings(self):
        summary = self.convert()
        folder = self.root / "run"
        folder.mkdir()
        (folder / "run_settings.json").write_text(json.dumps({"attempts": [
            {"n": 1, "started": "2026-10-22T18:24:05-04:00", "finished": "x", "versions": {"orthofinder": "3.1.5"},
             "coordinator_version": "0.13.0", "coordinator_version_note": "recorded afterwards",
             "settings": {"seed": 1}}]}))
        (folder / "tree_settings.json").write_text(json.dumps(
            {"model": "MFP", "since": "2026-10-30T10:00:00-04:00",
             "history": [{"model": "MFP;1000:TEST", "since": "2026-10-22T18:26:01-04:00"}]}))
        rec = C.run_record(folder / "run_settings.json", summary, "0.14.1", "0" * 32, "f" * 32, SPECIES)
        p = rec["parameters"]
        self.assertEqual((p["model_now"], p["model_since"]), ("MFP", "2026-10-30T10:00:00-04:00"))
        self.assertEqual(p["model_history"], [{"model": "MFP;1000:TEST", "since": "2026-10-22T18:26:01-04:00"}])
        self.assertEqual(p["coordinator_versions"], [{"attempt": 1, "version": "0.13.0", "note": "recorded afterwards"}])
