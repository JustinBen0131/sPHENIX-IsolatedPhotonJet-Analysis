"""Connected-chain checks on the small p+p / Au+Au fixture trees.

These are mechanism tests: the fixtures are synthetic engineering trees with
the older stored thresholds, so ``--thresholds-from tree`` is used.  They prove
that the TreeToHists stage reproduces the independent reference reducers
(``photonjet.analysis.purity`` and ``photonjet.analysis.reduce``) bin by bin,
that scoring preserves the tree contract and joins by candidate identity, and
that the purity-correction driver and the overlay run end to end.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import uproot
import yaml

import make_histograms
import score_trees
from photon_selection import load_config
from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection, recoil_histogram
from photonjet.io.tree_validation import validate

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
PYTHON = sys.executable


def open_config() -> dict:
    config = load_config(ROOT / "config" / "nominal.yaml")
    # The reference reducers apply no vertex or centrality window; match them.
    for system in ("pp", "auau"):
        config["systems"][system]["event"]["abs_vertex_z_max_cm"] = 1.0e9
        config["systems"][system]["event"].pop("centrality_percent", None)
    return config


class HistogramsMatchReferenceReducers(unittest.TestCase):
    def check_system(self, system: str) -> None:
        fixture = FIXTURES / f"photonjet_trees_{system}.root"
        config = open_config()
        payload = make_histograms.build([(fixture, None)], system=system, config=config, thresholds_from="tree", truth_signal_only=False)
        counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
        reference = purity_counts([fixture], non_tight_definition="bounded", isolation_radius=0.3)
        for index, region in enumerate("ABCD"):
            self.assertAlmostEqual(counts[index].sum(), reference["weighted_counts"][region], places=12, msg=f"{system} region {region}")
        xj_edges = config["recoil"]["xj_edges"]
        for spectrum_index, region in enumerate(("A", "C")):
            hist = recoil_histogram([fixture], RecoilSelection(region=region, isolation_radius=0.3), edges=xj_edges)
            spectra = np.asarray(payload["inclusive_recoil_spectra_ptgamma_xj"])[spectrum_index].sum(axis=0)
            sumw2 = np.asarray(payload["inclusive_recoil_sumw2_ptgamma_xj"])[spectrum_index].sum(axis=0)
            np.testing.assert_allclose(spectra, hist["sumw"], atol=1e-12, err_msg=f"{system} region {region} sumw")
            np.testing.assert_allclose(sumw2, hist["sumw2"], atol=1e-12, err_msg=f"{system} region {region} sumw2")

    def test_pp(self) -> None:
        self.check_system("pp")

    def test_auau(self) -> None:
        self.check_system("auau")

    def test_missing_working_points_stop_the_stage(self) -> None:
        with self.assertRaises(SystemExit):
            make_histograms.build([(FIXTURES / "photonjet_trees_pp.root", None)], system="pp", config=open_config(),
                                  thresholds_from="config", truth_signal_only=False)


class ScoringPreservesContract(unittest.TestCase):
    def test_flags_regenerate_from_config_and_join_by_identity(self) -> None:
        fixture = FIXTURES / "photonjet_trees_pp.root"
        config = open_config()
        config["systems"]["pp"]["working_points"].update({"tight_threshold": 0.5, "nontight_low": 0.1, "nontight_high": 0.4,
                                                          "isolated_max_gev": 4.0, "nonisolated_min_gev": 7.0})
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / fixture.name
            receipt = score_trees.score_file(fixture, output, system="pp", config=config, evaluator=None, model_path=None)
            self.assertEqual(receipt["score_source"], "stored")
            self.assertEqual(receipt["working_points_missing"], [])
            lines = validate([output], model_input_count=11)
            self.assertIn("tree_schema_exact=PASS", lines)
            self.assertIn("eventTree_flat_equivalence=PASS", lines)
            with uproot.open(output) as root:
                photons = root["photons"].arrays(["candidate_id_hi", "candidate_id_lo", "bdt_score", "bdt_tight_threshold", "bdt_is_tight", "iso_r03", "iso_r03_pass"], library="np")
                pairs = root["photonJets"].arrays(["candidate_id_hi", "candidate_id_lo", "bdt_tight_threshold", "bdt_is_tight", "iso_r03_pass"], library="np")
            np.testing.assert_allclose(photons["bdt_tight_threshold"], 0.5)
            np.testing.assert_array_equal(photons["bdt_is_tight"], (photons["bdt_score"] > 0.5).astype(np.int32))
            np.testing.assert_array_equal(photons["iso_r03_pass"], (photons["iso_r03"] < 4.0).astype(np.int32))
            by_candidate = {(int(h), int(l)): (int(t), int(i)) for h, l, t, i in zip(photons["candidate_id_hi"], photons["candidate_id_lo"], photons["bdt_is_tight"], photons["iso_r03_pass"])}
            for h, l, t, i in zip(pairs["candidate_id_hi"], pairs["candidate_id_lo"], pairs["bdt_is_tight"], pairs["iso_r03_pass"]):
                self.assertEqual(by_candidate[(int(h), int(l))], (int(t), int(i)))

    def test_duplicate_candidate_identity_is_rejected(self) -> None:
        fixture = FIXTURES / "photonjet_trees_pp.root"
        with tempfile.TemporaryDirectory() as temporary:
            with uproot.open(fixture) as root:
                trees = {name: root[name].arrays(library="np") for name in ("events", "photons", "jets", "photonJets", "truthPhotons", "truthJets", "recoTruthLinks")}
                event_tree = root["eventTree"].arrays(library="ak")
            trees["photons"]["candidate_id_lo"][1] = trees["photons"]["candidate_id_lo"][0]
            trees["photons"]["candidate_id_hi"][1] = trees["photons"]["candidate_id_hi"][0]
            tampered = Path(temporary) / "tampered.root"
            with uproot.recreate(tampered) as root:
                for name in ("events", "photons", "jets", "photonJets", "truthPhotons", "truthJets", "recoTruthLinks"):
                    root[name] = trees[name]
                root["eventTree"] = {field: event_tree[field] for field in event_tree.fields}
            with self.assertRaisesRegex(ValueError, "not unique"):
                score_trees.score_file(tampered, Path(temporary) / "out.root", system="pp", config=open_config(), evaluator=None, model_path=None)


class CorrectionDriverRuns(unittest.TestCase):
    def test_purity_correction_and_overlay(self) -> None:
        fixture = FIXTURES / "photonjet_trees_pp.root"
        config = open_config()
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            data = make_histograms.build([(fixture, None)], system="pp", config=config, thresholds_from="tree", truth_signal_only=False)
            # Synthetic leakage with the data's structure: region A dominated,
            # small B/C/D leakage.  The fixtures carry no nominal truth-signal
            # links, so this only exercises the driver and the kernels.
            leakage = copy.deepcopy(data)
            leakage["truth_signal_only"] = True
            counts = np.asarray(data["abcd_event_leading_counts_ptgamma"])
            leak = np.zeros_like(counts)
            leak[0] = np.maximum(counts[0], 1.0) * 10.0
            leak[1], leak[2], leak[3] = leak[0] * 0.05, leak[0] * 0.02, leak[0] * 0.001
            leakage["abcd_event_leading_counts_ptgamma"] = leak.tolist()
            leakage["abcd_event_leading_sumw2_ptgamma"] = leak.tolist()
            (work / "data.json").write_text(json.dumps(data), encoding="utf-8")
            (work / "leakage.json").write_text(json.dumps(leakage), encoding="utf-8")
            (work / "nominal.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
            completed = subprocess.run(
                [PYTHON, str(ROOT / "FinalAnalysis" / "run_corrections.py"), "--config", str(work / "nominal.yaml"), "--system", "pp",
                 "--data", str(work / "data.json"), "--leakage", str(work / "leakage.json"), "--output", str(work / "result.json")],
                capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads((work / "result.json").read_text(encoding="utf-8"))
            self.assertEqual(result["stage"], "purity_correction_only")
            self.assertEqual(len(result["purity"]), 3)
            self.assertEqual(np.asarray(result["density_per_ptgamma"]).shape, (3, len(config["recoil"]["xj_edges"]) - 1))
            overlay = subprocess.run(
                [PYTHON, str(ROOT / "FinalAnalysis" / "plot_overlay.py"), "--result", f"fixture={work / 'result.json'}", "--pt-bin", "0",
                 "--reference", f"ATLAS p+p (63-80 GeV)={ROOT / 'FinalAnalysis/reference/atlas_plb789_167_table1_xjgamma.csv'}:pp",
                 "--output", str(work / "overlay.png")],
                capture_output=True, text=True, check=False)
            self.assertEqual(overlay.returncode, 0, overlay.stderr)
            self.assertTrue((work / "overlay.png").is_file())


if __name__ == "__main__":
    unittest.main()
