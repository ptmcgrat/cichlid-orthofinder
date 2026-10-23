import gzip
import shutil
import tempfile
import unittest
from pathlib import Path

from cichlid_orthofinder import native

FIXTURE = Path(__file__).with_name("fixtures") / "pilot2" / "Results_Oct20"


class TestKeepList(unittest.TestCase):
    def test_keeps_the_results_and_leaves_out_the_search(self):
        with tempfile.TemporaryDirectory() as d:
            results = Path(d) / "Results_Oct20"
            shutil.copytree(FIXTURE, results)
            w = results / "WorkingDirectory"
            (w / "Species0.fa").write_text(">0_0\nMK\n")
            (w / "diamondDBSpecies0.dmnd").write_bytes(b"db")
            with gzip.open(w / "Blast0_0.txt.gz", "wb") as fh:
                fh.write(b"0_0\t0_0\t100\n")
            (w / "pickle").mkdir()
            (w / "pickle" / "x.pic").write_bytes(b"p")
            (w / "N0.ids.tsv").write_text("HOG\n")
            (w / "Sequences_ids").mkdir()
            sizes = native.curate(results, Path(d) / "native" / "our_species_tree")
            dest = Path(d) / "native" / "our_species_tree"
            kept = {p.relative_to(dest).as_posix() for p in dest.rglob("*")}
            for name in ("Orthogroups", "MultipleSequenceAlignments", "Resolved_Gene_Trees", "Log.txt",
                         "WorkingDirectory/SpeciesIDs.txt", "WorkingDirectory/SequenceIDs.txt",
                         "WorkingDirectory/N0.ids.tsv", "WorkingDirectory/clusters_OrthoFinder_I1.2.txt",
                         "WorkingDirectory/Trees_ids", "WorkingDirectory/Alignments_ids"):
                self.assertIn(name, kept)
            for name in ("WorkingDirectory/Species0.fa", "WorkingDirectory/diamondDBSpecies0.dmnd",
                         "WorkingDirectory/Blast0_0.txt.gz", "WorkingDirectory/pickle",
                         "WorkingDirectory/Sequences_ids",
                         "WorkingDirectory/clusters_OrthoFinder_I1.2.txt_id_pairs.txt"):
                self.assertNotIn(name, kept)
            self.assertEqual(set(sizes), {e.name for e in dest.iterdir()})
            self.assertTrue(all(n > 0 for n in sizes.values()))

    def test_refuses_an_existing_destination(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "dest").mkdir()
            with self.assertRaises(FileExistsError):
                native.curate(FIXTURE, Path(d) / "dest")


if __name__ == "__main__":
    unittest.main()
