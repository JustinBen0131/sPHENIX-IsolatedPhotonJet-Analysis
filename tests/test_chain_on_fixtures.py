"""Connected-chain checks on the small p+p / Au+Au fixture trees.

These are mechanism tests: the fixtures are synthetic engineering trees with
the older stored thresholds, so ``--thresholds-from tree`` is used.  They prove
that the TreeToHists stage reproduces the independent reference reducers
(``photonjet.analysis.purity`` and ``photonjet.analysis.reduce``) bin by bin,
that scoring preserves the tree contract and joins by candidate identity, and
that the correction/unfolding driver runs end to end on a real response bundle.
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

import make_histograms
import score_trees
from photon_selection import load_config
from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection, recoil_histogram
from photonjet.analysis.response_builder import ResponseBuildConfig, build_response
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
        payload = make_histograms.build([fixture], system=system, config=config, thresholds_from="tree", truth_signal_only=False)
        counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
        reference = purity_counts([fixture], non_tight_definition="bounded", isolation_radius=0.3)
        for index, region in enumerate("ABCD"):
            self.assertAlmostEqual(counts[index].sum(), reference["weighted_counts"][region], places=12, msg=f"{system} region {region}")
        xj_edges = config["recoil"]["xj_edges"]
        for spectrum_index, region in enumerate(("A", "C")):
            selection = RecoilSelection(region=region, isolation_radius=0.3)
            hist = recoil_histogram([fixture], selection, edges=xj_edges)
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
            make_histograms.build([FIXTURES / "photonjet_trees_pp.root"], system="pp", config=open_config(),
                                  thresholds_from="config", truth_signal_only=False)


class ScoringPreservesContract(unittest.TestCase):
    def test_flags_regenerate_from_config_and_join_by_identity(self) -> None:
        fixture = FIXTURES / "photonjet_trees_pp.root"
        config = open_config()
        wp = config["systems"]["pp"]["working_points"]
        wp.update({"tight_threshold": 0.5, "nontight_low": 0.1, "nontight_high": 0.4,
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


class CorrectionChainRuns(unittest.TestCase):
    """Purity correction on the fixtures, then unfolding on the same response."""

    def build_inputs(self, work: Path) -> tuple[Path, Path, dict]:
        fixture = FIXTURES / "photonjet_trees_pp.root"
        config = open_config()
        data = make_histograms.build([fixture], system="pp", config=config,
                                     thresholds_from="tree", truth_signal_only=False)
        # Prompt-photon leakage payload: region A dominated with small B/C/D
        # leakage.  Synthetic, because the fixtures carry no truth-signal links;
        # it exercises the ABCD solve and the correction arithmetic.
        leakage = copy.deepcopy(data)
        leakage["truth_signal_only"] = True
        counts = np.asarray(data["abcd_event_leading_counts_ptgamma"])
        leak = np.zeros_like(counts)
        leak[0] = np.maximum(counts[0], 1.0) * 10.0
        leak[1], leak[2], leak[3] = leak[0] * 0.05, leak[0] * 0.02, leak[0] * 0.001
        leakage["abcd_event_leading_counts_ptgamma"] = leak.tolist()
        leakage["abcd_event_leading_sumw2_ptgamma"] = leak.tolist()
        data_path, leakage_path = work / "data.json", work / "leakage.json"
        data_path.write_text(json.dumps(data), encoding="utf-8")
        leakage_path.write_text(json.dumps(leakage), encoding="utf-8")
        return data_path, leakage_path, config

    def test_purity_correction_runs_and_plots(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            data_path, leakage_path, config = self.build_inputs(work)
            config_path = work / "nominal.yaml"
            import yaml
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            completed = subprocess.run(
                [PYTHON, str(ROOT / "FinalAnalysis" / "run_corrections.py"),
                 "--data", str(data_path), "--leakage", str(leakage_path),
                 "--config", str(config_path), "--output", str(work / "result.json")],
                capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads((work / "result.json").read_text(encoding="utf-8"))
            self.assertEqual(len(result["abcd_solutions"]), 3)
            n_xj = len(config["recoil"]["xj_edges"]) - 1
            self.assertEqual(np.asarray(result["corrected_ptgamma_xj"]).shape, (3, n_xj))
            for solution in result["abcd_solutions"]:
                self.assertTrue(np.isfinite(solution["purity"]))
            overlay = subprocess.run(
                [PYTHON, str(ROOT / "FinalAnalysis" / "plot_overlay.py"),
                 "--result", f"fixture={work / 'result.json'}", "--stage", "corrected",
                 "--pt-bin", "0", "--output", str(work / "overlay.png")],
                capture_output=True, text=True, check=False)
            self.assertEqual(overlay.returncode, 0, overlay.stderr)
            self.assertTrue((work / "overlay.png").is_file())

    def test_unfolding_places_the_measured_spectrum_on_the_response_grid(self) -> None:
        """unfold() places the measured (pT, xJ) spectrum on the response grid.

        The response here is a synthetic diagonal one on the same axes the
        builder uses, because the two-event fixtures contain too few matched
        pairs to populate a response matrix.  The response arithmetic itself is
        covered by tests/test_response_unfolding.py and test_response_builder.py.
        """

        import run_corrections
        from photonjet.analysis.response import (
            CLASSIFICATION_RECO_PTGAMMA_EDGES,
            CLASSIFICATION_TRUTH_PTGAMMA_EDGES,
            COMMON_XJ_EDGES,
        )

        reco_edges = CLASSIFICATION_RECO_PTGAMMA_EDGES
        truth_edges = CLASSIFICATION_TRUTH_PTGAMMA_EDGES
        xj_edges = COMMON_XJ_EDGES
        n_xj = len(xj_edges) - 1
        n_truth, n_reco = (len(truth_edges) - 1) * n_xj, (len(reco_edges) - 1) * n_xj
        matrix = np.zeros((n_truth, n_reco))
        # Diagonal response: a truth (pT, xJ) cell migrates to the reconstructed
        # cell with the same pT low edge and the same xJ bin.
        for truth_pt in range(len(truth_edges) - 1):
            same = np.flatnonzero(np.isclose(reco_edges, truth_edges[truth_pt]))
            if not len(same) or same[0] >= len(reco_edges) - 1:
                continue
            reco_pt = int(same[0])
            for xj in range(n_xj):
                matrix[truth_pt * n_xj + xj, reco_pt * n_xj + xj] = 100.0
        response = {
            "reco_ptgamma_edges": reco_edges,
            "truth_ptgamma_edges": truth_edges,
            "xj_edges": xj_edges,
            "matrix": matrix,
            "misses": np.full(n_truth, 10.0),
            "unfolding_fakes": np.zeros(n_reco),
        }
        corrected = np.zeros((3, n_xj))
        corrected[:, n_xj // 2] = [10.0, 6.0, 3.0]
        unfolded = run_corrections.unfold(corrected, corrected.copy(), np.asarray(xj_edges),
                                          response, iterations=2, toys=3, seed=1)
        values = np.asarray(unfolded["unfolded_ptgamma_xj"])
        self.assertEqual(values.shape, (3, n_xj))
        self.assertTrue(np.all(np.isfinite(values)))
        # A diagonal response with a 100/(100+10) efficiency returns the
        # measured content scaled by 1/efficiency in the same xJ bin.
        np.testing.assert_allclose(values[:, n_xj // 2], np.asarray([10.0, 6.0, 3.0]) * 1.1, rtol=1e-9)
        self.assertTrue(np.isfinite(unfolded["refold_chi2_ndf"]))


if __name__ == "__main__":
    unittest.main()
