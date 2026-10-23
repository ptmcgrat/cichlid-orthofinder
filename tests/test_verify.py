import unittest

from cichlid_orthofinder import verify


class TestVerify(unittest.TestCase):
    def test_model_gate(self):
        g = verify.gate_model("full")
        self.assertTrue(g.ok)
        self.assertEqual(g.lines, ["model setting: MFP (full)"])
        g = verify.gate_model("hybrid")
        self.assertTrue(g.ok)
        self.assertIn("departs from the sheet", g.lines[0])
        self.assertFalse(verify.gate_model("MFP;1000").ok)

    def test_small_alignment_is_repeatable_and_aligned(self):
        a = verify.small_alignment()
        self.assertEqual(a, verify.small_alignment())
        rows = [l for l in a.splitlines() if not l.startswith(">")]
        self.assertEqual(len(rows), 6)
        self.assertEqual({len(r) for r in rows}, {120})

    def test_gate_printout(self):
        g = verify.Gate(3, "something")
        g.note("detail")
        out = []
        g.show(out.append)
        self.assertEqual(out, ["gate 3 FAIL  something", "         detail"])


class TestSpeciesGates(unittest.TestCase):
    """The species gates, on a stand-in for the coordinator's catalog (the real one is the gates' own run)."""

    def fake(self, d, files=("Aaa", "Bbb"), part_list=None, toml_list=None, tree_md5=None):
        import hashlib
        from pathlib import Path
        from types import SimpleNamespace
        d = Path(d)
        (d / "config").mkdir()
        tree = d / "config" / "tree.nwk"
        tree.write_text("(Aaa,Bbb);\n")
        species = list(files)
        (d / "config" / "releases.toml").write_text(
            'releases.v1.parts.OrthoFinder.species = [' + ", ".join(f'"{p}"' for p in (toml_list or species)) + ']\n')
        organisms = {}
        for p in species:
            f = d / f"{p}.faa"
            f.write_text(f">{p}_1\nMK\n")
            organisms[p] = SimpleNamespace(prefix=p, proteome_md5=hashlib.md5(f.read_bytes()).hexdigest())

        class Cat:
            meta = {"config": {"tree.nwk": tree_md5 or hashlib.md5(tree.read_bytes()).hexdigest()}}
            def species(self, release=None, part=None): return [organisms[p] for p in species]
            def part(self, release, part): return {"species": part_list or species, "species_tree": "tree.nwk"}
            def species_tree(self, release, part): return tree
            def __getitem__(self, p): return organisms[p]

        layout = SimpleNamespace(proteome=lambda p: SimpleNamespace(local=d / f"{p}.faa"))
        return Cat(), layout, tree

    def patched(self, cat, layout):
        from unittest import mock
        return mock.patch("cichlid_data.Catalog.load", return_value=cat), mock.patch("cichlid_data.get_layout", return_value=layout)

    def test_proteome_files_rehashed(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as d:
            cat, layout, _ = self.fake(d)
            a, b = self.patched(cat, layout)
            with a, b:
                self.assertTrue(verify.gate_proteome_files().ok)
                (Path(d) / "Bbb.faa").write_text(">Bbb_1\nMKX\n")          # a changed file
                g = verify.gate_proteome_files()
            self.assertFalse(g.ok)
            self.assertIn("Bbb: file", g.lines[1])

    def test_species_lists_print_the_difference(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            cat, layout, _ = self.fake(d, toml_list=["Aaa", "Ccc"])
            a, b = self.patched(cat, layout)
            with a, b:
                g = verify.gate_species_lists()
            self.assertFalse(g.ok)                                     # also not 35 species
            self.assertTrue(any("only in the catalog ['Bbb']" in l and "['Ccc']" in l for l in g.lines))

    def test_species_tree_checksum(self):
        import json
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as d:
            cat, layout, tree = self.fake(d)
            run = Path(d) / "work" / "v1" / "our_species_tree"
            run.mkdir(parents=True)
            (run / "run_settings.json").write_text(json.dumps(
                {"species_tree_text": "(Aaa,Bbb);", "attempts": [{"settings": {"species_tree": str(tree)}}]}))
            a, b = self.patched(cat, layout)
            with a, b:
                self.assertTrue(verify.gate_species_tree(Path(d) / "work").ok)
                tree.write_text("(Bbb,Aaa);\n")
                self.assertFalse(verify.gate_species_tree(Path(d) / "work").ok)
