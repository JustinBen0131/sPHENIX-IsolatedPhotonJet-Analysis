"""Frontend contracts only: no ROOT, Fun4All, production inputs or science tests."""

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch


TREE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("frontend", TREE / "scripts" / "produce.py")
frontend = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frontend)

COMMANDS = {
    "runPP": ["pp_data"],
    "runAuAu": ["auau_data"],
    "runPhotonJetSim": ["pp_photon_sim"],
    "runInclusiveJetSim": ["pp_inclusive_sim"],
    "runEmbeddedPhotonJetSim": ["auau_photon_embedded"],
    "runEmbeddedInclusiveJetSim": ["auau_inclusive_embedded"],
}
COMMANDS["runAllData"] = COMMANDS["runPP"] + COMMANDS["runAuAu"]
COMMANDS["runAllSim"] = sum(list(COMMANDS.values())[2:6], [])
COMMANDS["runAll"] = COMMANDS["runAllData"] + COMMANDS["runAllSim"]


class FrontendTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="tree frontend tests ")
        self.addCleanup(self.temporary.cleanup)
        self.tree = Path(self.temporary.name) / "Tree Production"
        for directory in ("config", "scripts", "macros", "fixtures"):
            (self.tree / directory).mkdir(parents=True, exist_ok=True)
        for name in ("produce_trees.sh", "scripts/produce.py", "config/production_populations.json",
                     "config/tree_production.yaml", "macros/Fun4All_PhotonJetTree.C"):
            shutil.copyfile(TREE / name, self.tree / name)
        self.binding_path = self.tree / "config" / "production_populations.json"
        self.settings = json.loads(self.binding_path.read_text())
        self.settings["campaign"] = "unit-fixture-not-production"
        self.manifest = self.tree / "fixtures" / "population.json"
        inputs = self.tree / "fixtures" / "source.list"
        inputs.write_text("METADATA-FIXTURE-NEVER-EXECUTED\n")
        self.source = dict(source_file_ordinal=7, run=1234, segment=0,
                           first_entry=0, event_count=20, sample="unit fixture", period="", si_di_role="",
                           input_list="source.list", input_list_sha256=frontend.sha256(inputs))
        self.population = dict(schema_version=1, campaign=self.settings["campaign"],
                               lane="pp_data", sources=[self.source])
        self.args = SimpleNamespace(events=None, source=None, dry_run=True)
        for name, value in (("TREE_DIR", self.tree), ("BINDINGS", self.binding_path),
                            ("MACRO", self.tree / "macros" / "Fun4All_PhotonJetTree.C")):
            manager = patch.object(frontend, name, value)
            manager.start()
            self.addCleanup(manager.stop)

    def bind(self):
        self.manifest.write_text(json.dumps(self.population))
        self.settings["lanes"]["pp_data"].update(
            manifest="../fixtures/population.json", manifest_sha256=frontend.sha256(self.manifest))
        self.binding_path.write_text(json.dumps(self.settings))

    def plan(self):
        self.bind()
        return frontend.plan_population("pp_data", self.settings["lanes"]["pp_data"], self.settings, self.args)

    def cli(self, *arguments):
        return subprocess.run(["/bin/bash", str(self.tree / "produce_trees.sh"), *arguments],
                              cwd=self.temporary.name, text=True, capture_output=True)

    def test_help_and_no_arguments_match_from_another_directory(self):
        help_result = self.cli("--help")
        self.assertEqual(help_result.returncode, 0)
        self.assertEqual(self.cli().stdout, help_result.stdout)
        for command in COMMANDS:
            self.assertIn(command, help_result.stdout)

    def test_all_nine_commands_fail_closed_and_aggregates_keep_order(self):
        for command, lanes in COMMANDS.items():
            with self.subTest(command=command):
                result = self.cli(command)
                self.assertNotEqual(result.returncode, 0)
                positions = [result.stderr.index(f"lanes.{lane}.manifest") for lane in lanes]
                self.assertEqual(positions, sorted(positions))
                self.assertEqual(result.stderr.count("not yet configured"), len(lanes))
        self.assertFalse((self.tree / "output").exists())

    def test_invalid_commands_options_and_event_ranges(self):
        for arguments in (("runDI",), ("runPP", "--submit"), ("runPP", "--events"),
                          ("runPP", "--events", "0"), ("runPP", "--events", "-1"),
                          ("runPP", "--events", "2147483648"), ("runAll", "--source", "7")):
            with self.subTest(arguments=arguments):
                self.assertNotEqual(self.cli(*arguments).returncode, 0)
        self.assertIn("Usage:", self.cli("unknown").stderr)

    def test_dry_run_does_not_invoke_root_or_create_output(self):
        self.bind()
        with patch.object(frontend.subprocess, "run") as run:
            with patch.object(frontend.sys, "argv", ["produce.py", "--dry-run", "pp_data"]):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    self.assertEqual(frontend.main(), 0)
            run.assert_not_called()
        self.assertIn("DRY RUN", output.getvalue())
        self.assertFalse((self.tree / "output").exists())
        self.assertEqual(self.cli("runPP", "--dry-run").returncode, 0)

    def test_aggregate_preflights_missing_lane_before_any_execution(self):
        self.bind()
        with patch.object(frontend, "execute_source") as execute:
            with patch.object(frontend.sys, "argv", ["produce.py", "pp_data", "auau_data"]):
                with self.assertRaisesRegex(frontend.FrontendError, "auau_data.*not yet configured"):
                    frontend.main()
            execute.assert_not_called()
        self.assertFalse((self.tree / "output").exists())

    def test_development_preserves_identity_and_uses_distinct_output(self):
        nominal = self.plan()
        self.args.events, self.args.source = 10, 7
        restricted = self.plan()
        self.assertEqual(restricted["jobs"][0]["events"], 10)
        self.assertEqual(restricted["jobs"][0]["source_manifest_sha256"], nominal["manifest_sha256"])
        self.assertEqual(restricted["jobs"][0]["source_file_ordinal"], 7)
        self.assertNotEqual(restricted["output"], nominal["output"])
        self.args.events = 30
        self.assertEqual(self.plan()["jobs"][0]["events"], 20)
        self.args.source = 0
        with self.assertRaisesRegex(frontend.FrontendError, "does not exist"):
            self.plan()

    def test_empty_duplicate_missing_and_unknown_metadata_are_rejected(self):
        variants = [[], [dict(self.source, run=None)], [self.source, self.source],
                    [dict(self.source, event_count=True)], [dict(self.source, unexpected="ignored?")],
                    [dict(self.source, production_config_sha256="a" * 64)]]
        for sources in variants:
            with self.subTest(sources=sources):
                self.population["sources"] = sources
                with self.assertRaises(frontend.FrontendError):
                    self.plan()

    def test_digest_and_lane_mismatch_are_rejected(self):
        self.bind()
        self.manifest.write_text(self.manifest.read_text() + "\n")
        with self.assertRaisesRegex(frontend.FrontendError, "SHA256 mismatch"):
            frontend.plan_population("pp_data", self.settings["lanes"]["pp_data"], self.settings, self.args)
        self.population["lane"] = "auau_data"
        with self.assertRaisesRegex(frontend.FrontendError, "lane mismatch"):
            self.plan()

    def test_existing_output_refused(self):
        plan = self.plan()
        Path(plan["output"]).mkdir(parents=True)
        with self.assertRaisesRegex(frontend.FrontendError, "refusing overwrite"):
            self.plan()

    def test_macro_arguments_are_quoted_without_a_shell(self):
        self.source["sample"] = 'quote" and backslash\\ and $(do-not-run)'
        job = self.plan()["jobs"][0]
        invocation = frontend.macro_argument(job)
        self.assertIn(json.dumps(self.source["sample"]), invocation)
        self.assertTrue(invocation.startswith(str(frontend.MACRO) + "("))
        self.assertIn(',1234,0,20,0,', invocation)

    def test_child_failure_propagates_and_never_checks_or_certifies_output(self):
        job = self.plan()["jobs"][0]
        Path(job["output_file"]).parent.parent.mkdir(parents=True)
        with patch.object(frontend.subprocess, "run", return_value=SimpleNamespace(returncode=17)) as run:
            with patch.object(frontend, "check_completion") as check:
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    status = frontend.execute_source(job, "root", {}, None, "a" * 64)
                self.assertEqual(status, 17)
                check.assert_not_called()
                self.assertNotIn("shell", run.call_args.kwargs)
                self.assertEqual(run.call_args.args[0][:4], ["root", "-l", "-b", "-q"])
        self.assertFalse(Path(job["output_file"]).with_name("completion.json").exists())

    def test_zero_exit_still_requires_terminal_metadata(self):
        job = self.plan()["jobs"][0]
        Path(job["output_file"]).parent.parent.mkdir(parents=True)
        with patch.object(frontend.subprocess, "run", return_value=SimpleNamespace(returncode=0)):
            with patch.object(frontend, "check_completion", side_effect=ValueError("missing completion")):
                with contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaisesRegex(frontend.FrontendError, "completion check failed"):
                        frontend.execute_source(job, "root", {}, None, "a" * 64)
        self.assertFalse(Path(job["output_file"]).with_name("completion.json").exists())

    def product(self, job):
        metadata = {key: str(job[key]) for key in ("run", "segment", "source_file_ordinal", "first_entry",
                    "source_manifest_sha256", "sample", "period", "si_di_role", "configuration_sha256")}
        metadata.update(max_events=str(job["events"]), macro_sha256=frontend.sha256(frontend.MACRO),
                        producer_library_sha256="a" * 64)
        return {"metadata": "\n".join(f"{k}={v}" for k, v in metadata.items()),
                "completion": "completion_status=complete\nupstream_observed_events=20\n"
                              "encountered_events=19\nretained_events=19\nupstream_rejected_events=1",
                "Events": SimpleNamespace(num_entries=19),
                "UpstreamRejectedEvents": SimpleNamespace(num_entries=1),
                "Sources": Mock(num_entries=1)}

    def test_completion_rejects_aborted_short_or_mismatched_products(self):
        job = self.plan()["jobs"][0]
        product = self.product(job)
        # Mock only the reader boundary; this does not construct a ROOT file.
        source = Mock()
        source.array.return_value = [True]
        class SourceTable(dict):
            num_entries = 1
        product["Sources"] = SourceTable(completed=source)
        reader = Mock()
        reader.open.side_effect = lambda _: contextlib.nullcontext(product)
        frontend.check_completion(job, reader, "a" * 64)
        good = product["completion"]
        for bad in (good.replace("complete", "aborted"), good.replace("observed_events=20", "observed_events=19")):
            product["completion"] = bad
            with self.assertRaises(frontend.FrontendError):
                frontend.check_completion(job, reader, "a" * 64)
        product["completion"] = good
        with self.assertRaisesRegex(frontend.FrontendError, "provenance mismatch"):
            frontend.check_completion(job, reader, "b" * 64)


if __name__ == "__main__":
    unittest.main()
