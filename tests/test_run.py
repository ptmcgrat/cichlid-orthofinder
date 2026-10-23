import subprocess
import tempfile
import unittest
from pathlib import Path

from cichlid_orthofinder import run


class TestBuildCommand(unittest.TestCase):
    def test_with_species_tree(self):
        cmd = run.build_command("in", "out", "tree.nwk", threads=32)
        self.assertEqual(cmd, ["cichlid-orthofinder-of", "-f", "in", "-o", "out", "-M", "msa", "-T", "cichlid_tree",
                               "-s", "tree.nwk", "-t", "32"])

    def test_own_species_tree_has_no_s(self):
        cmd = run.build_command("in", "out", None)
        self.assertNotIn("-s", cmd)
        self.assertEqual(cmd[-4:], ["-M", "msa", "-T", "cichlid_tree"])

    def test_resume_starts_from_the_search(self):
        cmd = run.resume_command("out/Results_Oct22/WorkingDirectory", "tree.nwk", threads=8)
        self.assertEqual(cmd[:3], ["cichlid-orthofinder-of", "-b", "out/Results_Oct22/WorkingDirectory"])
        self.assertIn("cichlid_tree", cmd)
        self.assertNotIn("-f", cmd)


class TestRegister(unittest.TestCase):
    def test_adds_entry_and_keeps_others(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / "config_orthofinder_user.json"
            cfg.write_text(json.dumps({"iqtree3_mt": {"program_type": "tree"}}))
            run.register_tree_program(cfg)
            data = json.loads(cfg.read_text())
            self.assertEqual(data["cichlid_tree"], run.TREE_ENTRY)
            self.assertIn("iqtree3_mt", data)
            before = cfg.stat().st_mtime_ns
            run.register_tree_program(cfg)                      # already right: not rewritten
            self.assertEqual(cfg.stat().st_mtime_ns, before)


class TestLatestResults(unittest.TestCase):
    def test_follows_resumes(self):
        with tempfile.TemporaryDirectory() as d:
            first = Path(d) / "Results_Oct22"
            (first / "WorkingDirectory").mkdir(parents=True)
            self.assertEqual(run.latest_results(d), first)
            second = first / "WorkingDirectory" / "OrthoFinder" / "Results_Oct22"
            (second / "WorkingDirectory").mkdir(parents=True)
            self.assertEqual(run.latest_results(d), second)


class TestPrepareInput(unittest.TestCase):
    def test_links_named_by_prefix_in_order(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            src = d / "src"
            src.mkdir()
            for p in ("Mzeb", "Acal"):
                (src / f"{p}.faa").write_text(f">{p}_1\nM\n")
            folder = run.prepare_input(d / "input", {"Mzeb": src / "Mzeb.faa", "Acal": src / "Acal.faa"})
            self.assertEqual(sorted(x.name for x in folder.iterdir()), ["Acal.faa", "Mzeb.faa"])
            self.assertTrue((folder / "Mzeb.faa").is_symlink())
            self.assertEqual((folder / "Mzeb.faa").resolve(), (src / "Mzeb.faa").resolve())

    def test_refuses_non_empty_folder_and_missing_proteome(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "input").mkdir()
            (d / "input" / "x").write_text("")
            with self.assertRaises(RuntimeError):
                run.prepare_input(d / "input", {})
            with self.assertRaises(FileNotFoundError):
                run.prepare_input(d / "new", {"Mzeb": d / "missing.faa"})


class TestFindResults(unittest.TestCase):
    def test_one_results_folder(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "Results_Oct20").mkdir()
            self.assertEqual(run.find_results(d).name, "Results_Oct20")
            (Path(d) / "Results_Oct21").mkdir()
            with self.assertRaises(RuntimeError):
                run.find_results(d)


class TestVersions(unittest.TestCase):
    def test_parses_each_program(self):
        outputs = {"orthofinder": "\x1b[1mOrthoFinder:v3.1.5\x1b[0m\n",
                   "iqtree3": "IQ-TREE version 3.1.3 for Linux x86 64-bit\n",
                   "famsa": "FAMSA\n  version 2.4.1-45c9b2b (2025-05-09)\n",
                   "diamond": "diamond version 2.2.8\n",
                   "mcl": "mcl 22-282\n"}

        def fake(cmd, **kw):
            return subprocess.CompletedProcess(cmd, 0, stdout=outputs[cmd[0]], stderr="")

        self.assertEqual(run.program_versions(fake), {"orthofinder": "3.1.5", "iqtree3": "3.1.3",
                                                      "famsa": "2.4.1-45c9b2b", "diamond": "2.2.8",
                                                      "mcl": "22-282"})

    def test_missing_program_is_none(self):
        def fake(cmd, **kw):
            raise FileNotFoundError(cmd[0])

        self.assertTrue(all(v is None for v in run.program_versions(fake).values()))


if __name__ == "__main__":
    unittest.main()


class TestModelSetting(unittest.TestCase):
    def test_shortcuts_and_history(self):
        import json, tempfile
        from cichlid_orthofinder import run as R
        self.assertEqual(R.model_rule("hybrid"), "MFP;1000:TEST")
        self.assertEqual(R.model_rule("FULL"), "MFP")
        self.assertEqual(R.model_rule(R.DEFAULT_MODEL), "MFP", "the default is the sheet's setting")
        with self.assertRaises(ValueError):
            R.model_rule("MFP;1000")
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(R.read_model(d))
            R.write_model(d, "hybrid")
            R.write_model(d, "full")
            self.assertEqual(R.read_model(d), "MFP")
            saved = json.loads((Path(d) / R.LIVE_SETTINGS).read_text())
            self.assertEqual([h["model"] for h in saved["history"]], ["MFP;1000:TEST"])

    def test_environment_points_trees_at_the_live_file(self):
        from cichlid_orthofinder import run as R
        env = R.tree_env("/c", 1, "1", 32, "ascending", "MFP;1000:TEST", "/w/tree_settings.json")
        self.assertEqual(env["CICHLID_TREE_MODEL"], "MFP;1000:TEST")
        self.assertEqual(env["CICHLID_TREE_SETTINGS"], "/w/tree_settings.json")

    def test_set_model_command(self):
        import io, tempfile
        from contextlib import redirect_stdout
        from cichlid_orthofinder import cli, run as R
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(cli.main(["set-model", d, "full"]), 1)      # not a run's folder
            (Path(d) / "run_settings.json").write_text("{}")
            out = io.StringIO()
            with redirect_stdout(out):
                self.assertEqual(cli.main(["set-model", d, "full"]), 0)
                self.assertEqual(cli.main(["set-model", d]), 0)
            self.assertEqual(R.read_model(d), "MFP")
            self.assertIn("model setting: MFP (full)", out.getvalue())
            self.assertEqual(cli.main(["set-model", d, "MFP;x:TEST"]), 1)
            self.assertEqual(R.read_model(d), "MFP")                      # a bad rule changes nothing


class TestSearchCheck(unittest.TestCase):
    """A run killed during the all-against-all search must not be resumed from partial hits."""

    def make(self, d, n=3, skip=(), cut=()):
        import gzip
        wd = Path(d) / "WorkingDirectory"
        wd.mkdir()
        (wd / "SpeciesIDs.txt").write_text("".join(f"{i}: S{i}.faa\n" for i in range(n)))
        for i in range(n):
            for j in range(n):
                f = wd / f"Blast{i}_{j}.txt.gz"
                if (i, j) in skip:
                    continue
                with gzip.open(f, "wb") as fh:
                    fh.write(b"0_1\t1_2\t90.0\n" * 2000)
                if (i, j) in cut:                                  # killed while writing
                    f.write_bytes(f.read_bytes()[: f.stat().st_size // 2])
        return wd

    def test_complete_search(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(run.search_problems(self.make(d)), [])

    def test_missing_and_cut_files_are_found(self):
        with tempfile.TemporaryDirectory() as d:
            problems = run.search_problems(self.make(d, skip={(2, 1)}, cut={(0, 2)}))
            self.assertEqual(len(problems), 2)
            self.assertTrue(any("Blast2_1.txt.gz missing" in p for p in problems))
            self.assertTrue(any("Blast0_2.txt.gz unreadable" in p for p in problems))

    def test_search_never_started(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(run.search_problems(Path(d)), ["SpeciesIDs.txt is missing (the search never started)"])

    def test_set_aside_keeps_the_folder(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "output"
            (out / "Results_X").mkdir(parents=True)
            aside = run.set_aside(out)
            self.assertFalse(out.exists())
            self.assertTrue((aside / "Results_X").is_dir())
            self.assertTrue(aside.name.startswith("output_incomplete_"))


class TestSearchReuse(unittest.TestCase):
    """A saved search is reused only if it is complete and of exactly the same proteomes."""

    def make(self, d):
        import gzip
        d = Path(d)
        prot = {"Aaa": d / "Aaa.faa", "Bbb": d / "Bbb.faa"}
        prot["Aaa"].write_text(">Aaa_1 desc\nMKV\nLL*\n>Aaa_2\nMRR\n")
        prot["Bbb"].write_text(">Bbb_7\nmqq\n")
        wd = d / "WorkingDirectory"
        wd.mkdir()
        (wd / "SpeciesIDs.txt").write_text("0: Aaa.faa\n1: Bbb.faa\n")
        (wd / "SequenceIDs.txt").write_text("0_0: Aaa_1 desc\n0_1: Aaa_2\n1_0: Bbb_7\n")
        (wd / "Species0.fa").write_text(">0_0\nMKVLL\n>0_1\nMRR\n")
        (wd / "Species1.fa").write_text(">1_0\nMQQ\n")
        for i in range(2):
            for j in range(2):
                with gzip.open(wd / f"Blast{i}_{j}.txt.gz", "wb") as fh:
                    fh.write(b"0_0\t1_0\t50\n")
        (wd / "diamondDBSpecies0.dmnd").write_bytes(b"db")
        (wd / "clusters_OrthoFinder_I1.5.txt").write_text("old clusters\n")
        (wd / "Trees_ids").mkdir()
        return wd, prot

    def test_same_proteomes_match(self):
        with tempfile.TemporaryDirectory() as d:
            wd, prot = self.make(d)
            self.assertEqual(run.search_mismatches(wd, prot), [])

    def test_changed_sequence_is_caught(self):
        with tempfile.TemporaryDirectory() as d:
            wd, prot = self.make(d)
            prot["Aaa"].write_text(">Aaa_1 desc\nMKV\nLA\n>Aaa_2\nMRR\n")
            self.assertIn("Aaa: sequence 1", run.search_mismatches(wd, prot)[0])

    def test_changed_id_order_count_and_species_are_caught(self):
        with tempfile.TemporaryDirectory() as d:
            wd, prot = self.make(d)
            prot["Aaa"].write_text(">Aaa_2\nMRR\n>Aaa_1\nMKVLL\n")
            self.assertTrue(run.search_mismatches(wd, prot))
            prot["Aaa"].write_text(">Aaa_1\nMKVLL\n")
            self.assertIn("1 in the proteome", run.search_mismatches(wd, prot)[0])
            self.assertIn("species differ", run.search_mismatches(wd, {"Aaa": prot["Aaa"]})[0])

    def test_copy_takes_only_the_search_and_leaves_the_original(self):
        with tempfile.TemporaryDirectory() as d:
            wd, prot = self.make(d)
            before = sorted(p.name for p in wd.iterdir())
            dest = run.copy_search(wd, Path(d) / "new" / "WorkingDirectory")
            self.assertEqual(sorted(p.name for p in dest.iterdir()),
                             ["Blast0_0.txt.gz", "Blast0_1.txt.gz", "Blast1_0.txt.gz", "Blast1_1.txt.gz",
                              "SequenceIDs.txt", "Species0.fa", "Species1.fa", "SpeciesIDs.txt",
                              "diamondDBSpecies0.dmnd"])
            self.assertEqual(sorted(p.name for p in wd.iterdir()), before)
            self.assertEqual(run.search_problems(dest), [])
            self.assertEqual(run.search_mismatches(dest, prot), [])


class TestSection2Gates(unittest.TestCase):
    """Brief section 2: a recorded model setting (default full); every program version; the code's commit."""

    def test_model_defaults_to_full_and_is_printed(self):
        from unittest import mock
        for model, expect, departs in ((None, "model setting: MFP (full)", False),
                                       ("hybrid", "model setting: MFP;1000:TEST (hybrid)", True)):
            said = []
            with tempfile.TemporaryDirectory() as d, \
                    mock.patch.object(run, "program_versions", return_value={"orthofinder": None}), \
                    mock.patch.object(run.subprocess, "run", return_value=subprocess.CompletedProcess([], 0)), \
                    mock.patch("cichlid_data.get_store"), mock.patch("cichlid_data.Catalog"), \
                    mock.patch("cichlid_data.get_layout"):
                with self.assertRaises(RuntimeError):            # stops at the versions, after the model line
                    run.run(d, run_name="own_species_tree", model=model, dry_run=True, say=said.append)
            line = next(l for l in said if l.startswith("model setting:"))
            self.assertTrue(line.startswith(expect), line)
            self.assertEqual("departs from the sheet" in line, departs, line)

    def test_refuses_a_missing_version(self):
        from unittest import mock
        with tempfile.TemporaryDirectory() as d, \
                mock.patch.object(run, "program_versions", return_value={"orthofinder": "3.1.5", "famsa": None}), \
                mock.patch.object(run.subprocess, "run", return_value=subprocess.CompletedProcess([], 0)), \
                mock.patch("cichlid_data.get_store"), mock.patch("cichlid_data.Catalog"), \
                mock.patch("cichlid_data.get_layout"):
            with self.assertRaisesRegex(RuntimeError, "could not read the version of famsa"):
                run.run(d, run_name="own_species_tree", model="hybrid", dry_run=True)

    def test_code_record(self):
        rec = run.code_record(config="/nonexistent/config.json")
        self.assertIn("run.py", rec["source_sha256"])
        self.assertEqual(len(rec["source_sha256"]["run.py"]), 64)
        self.assertEqual(rec["tree_program"], run.TREE_ENTRY)
        self.assertIsNone(rec["orthofinder_user_config_entry"])
        self.assertEqual(rec["launcher_patch_md5"], "e875a68a2dfbd6a5c390d5a28d6d8cd4")
        t = rec["iqtree_argument_template"]
        self.assertEqual(t[0], "iqtree3")
        self.assertIn("-seed", t)
        self.assertEqual(t[t.index("-m") + 1], "<model>")

    def test_cli_model_defaults_to_full(self):
        from unittest import mock
        from cichlid_orthofinder import cli
        with mock.patch.object(run, "run") as fake:
            self.assertEqual(cli.main(["run", "--work", "/nonexistent", "--dry-run"]), 0)
            self.assertEqual(fake.call_args.kwargs["model"], "full")
            cli.main(["run", "--work", "/nonexistent", "--model", "hybrid", "--dry-run"])
            self.assertEqual(fake.call_args.kwargs["model"], "hybrid")

    def test_missing_versions(self):
        self.assertEqual(run.missing_versions({"orthofinder": "3.1.5", "famsa": None, "mcl": ""}), ["famsa", "mcl"])
        self.assertEqual(run.missing_versions({"orthofinder": "3.1.5"}), [])

    def test_version_patterns_read_each_program(self):
        outputs = {"orthofinder": "OrthoFinder:v3.1.5\n", "iqtree3": "IQ-TREE version 3.1.3 for Linux",
                   "famsa": "FAMSA (Fast and Accurate Multiple Sequence Alignment)\n  version 2.4.1-45c9b2b (2025-05-09)",
                   "diamond": "diamond version 2.2.8", "mcl": "mcl 22-282"}

        def fake(cmd, **kw):
            return subprocess.CompletedProcess(cmd, 0, stdout=outputs[cmd[0]], stderr="")
        self.assertEqual(run.program_versions(fake), {"orthofinder": "3.1.5", "iqtree3": "3.1.3",
                                                      "famsa": "2.4.1-45c9b2b", "diamond": "2.2.8", "mcl": "22-282"})

    def test_repo_state(self):
        def fake(cmd, **kw):
            out = "abc123\n" if "rev-parse" in cmd else " M src/x.py\n"
            return subprocess.CompletedProcess(cmd, 0, stdout=out, stderr="")
        self.assertEqual(run.repo_state(fake), {"commit": "abc123", "dirty": True})

        def no_git(cmd, **kw):
            return subprocess.CompletedProcess(cmd, 128, stdout="", stderr="fatal")
        self.assertEqual(run.repo_state(no_git), {"commit": None, "dirty": None})

    def test_launcher_checksum_is_the_briefs(self):
        self.assertEqual(run.TREES_MSA_MD5, "e875a68a2dfbd6a5c390d5a28d6d8cd4")


class TestCoordinatorVersionPerAttempt(unittest.TestCase):
    def test_earlier_attempts_get_the_version_they_ran_under(self):
        attempts = [{"n": 1, "started": "2026-10-22T18:24:05-04:00"},
                    {"n": 2, "started": "2026-11-01T09:00:00-04:00"},
                    {"n": 3, "started": "2026-11-02T09:00:00-04:00", "coordinator_version": "0.14.1"}]
        self.assertEqual(run.note_earlier_coordinator(attempts), 2)
        self.assertEqual(attempts[0]["coordinator_version"], "0.13.0")
        self.assertIn("recorded afterwards", attempts[0]["coordinator_version_note"])
        self.assertIsNone(attempts[1]["coordinator_version"])
        self.assertEqual(attempts[1]["coordinator_version_note"], "not recorded by this attempt")
        self.assertEqual(attempts[2], {"n": 3, "started": "2026-11-02T09:00:00-04:00", "coordinator_version": "0.14.1"})
        self.assertEqual(run.note_earlier_coordinator(attempts), 0)         # only once


class TestCliPassthrough(unittest.TestCase):
    def test_convert_native_verify_reach_their_modules(self):
        from unittest import mock
        from cichlid_orthofinder import cli
        for name in ("convert", "native", "verify"):
            with mock.patch(f"cichlid_orthofinder.{name}.main", return_value=7) as m:
                self.assertEqual(cli.main([name, "--x", "1"]), 7)
                m.assert_called_once_with(["--x", "1"])
