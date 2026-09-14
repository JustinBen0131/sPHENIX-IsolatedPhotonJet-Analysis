"""Per-event photon response, leakage exclusivity and weighted histograms on synthetic trees."""

from __future__ import annotations

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

import make_histograms
from photon_selection import load_config, truth_signal_mask
from photonjet.analysis.photon_response import build_photon_response
from photonjet.analysis.reduce import RecoilSelection
from sample_weights import load_manifest
from synthetic_trees import match_link, signal_truth_photon, write_tree_file

ROOT = Path(__file__).resolve().parents[1]
CONFIG = load_config(ROOT / "config" / "nominal.yaml")
MANIFEST = load_manifest(ROOT / "config" / "samples.yaml")
SELECTION = RecoilSelection(region="A", isolation_radius=0.3)


def truth_signal(columns):
    return truth_signal_mask(columns, CONFIG["truth_signal"])


def event(number: int, **overrides):
    row = {"source_file_index": 0, "event_id_hi": 1, "event_id_lo": number, "run": 1, "event_weight": 1.0,
           "centrality": -1.0, "vertex_z": 5.0, "terminal_status": 0}
    row.update(overrides)
    return row


def photon(number: int, candidate: int, et: float, *, score: float = 0.9, iso: float = 1.0, ordinal: int = 0):
    return {"source_file_index": 0, "event_id_hi": 1, "event_id_lo": number, "candidate_id_hi": 5, "candidate_id_lo": candidate,
            "photon_encounter_ordinal": ordinal, "photon_et": et, "photon_eta": 0.1, "photon_phi": 0.2,
            "bdt_score": score, "bdt_tight_threshold": 0.5, "bdt_nontight_low_threshold": 0.1, "bdt_nontight_high_threshold": 0.4,
            "iso_r03": iso, "iso_r03_threshold": 4.0, "iso_r03_nonisolated_threshold": 7.0,
            "iso_r04": iso, "iso_r04_threshold": 4.0, "iso_r04_nonisolated_threshold": 7.0, "bdt_input_count": 11}


class PhotonResponseTest(unittest.TestCase):
    def test_fill_placement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = write_tree_file(
                Path(temporary) / "sim.root",
                events=[event(n) for n in range(1, 7)],
                photons=[photon(1, 11, 21.0), photon(2, 12, 12.0), photon(3, 13, 18.0), photon(4, 14, 16.0), photon(5, 15, 30.0)],
                truth_photons=[
                    signal_truth_photon(0, 1, 1, 9, 1, 22.0), signal_truth_photon(0, 1, 2, 9, 2, 30.0),
                    signal_truth_photon(0, 1, 4, 9, 4, 8.0), signal_truth_photon(0, 1, 5, 9, 5, 45.0),
                    signal_truth_photon(0, 1, 6, 9, 6, 20.0),
                ],
                links=[match_link(0, 1, 1, 5, 11, 9, 1), match_link(0, 1, 2, 5, 12, 9, 2),
                       match_link(0, 1, 4, 5, 14, 9, 4), match_link(0, 1, 5, 5, 15, 9, 5)],
            )
            response = build_photon_response([path], system="pp", selection=SELECTION, truth_signal=truth_signal)
        # truth edges 5,10,15,20,25,35,40 ; reco edges 10,15,20,25,35,40
        self.assertEqual(response.matrix[3, 2], 1.0)     # truth 22 -> reco 21
        self.assertEqual(response.matrix[4, 0], 1.0)     # truth 30 -> reco 12
        self.assertEqual(response.matrix[0, 1], 1.0)     # truth 8 -> reco 16
        self.assertEqual(response.fakes[1], 1.0)         # reco 18 without truth
        self.assertEqual(response.boundary_fakes[3], 1.0)  # truth 45 off grid, reco 30
        self.assertEqual(response.misses[3], 1.0)        # truth 20 without reco
        self.assertEqual(response.provenance["counts"]["matched"], 4)   # events 1, 2, 4 and the off-grid 5
        response.check_conservation()

    def test_two_signal_truth_photons_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = write_tree_file(Path(temporary) / "sim.root", events=[event(1)], photons=[photon(1, 11, 21.0)],
                                   truth_photons=[signal_truth_photon(0, 1, 1, 9, 1, 22.0), signal_truth_photon(0, 1, 1, 9, 2, 25.0)])
            with self.assertRaisesRegex(ValueError, "more than one truth-signal photon"):
                build_photon_response([path], system="pp", selection=SELECTION, truth_signal=truth_signal)

    def test_weight_override_and_contract_rejections(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = write_tree_file(
                Path(temporary) / "sim.root", events=[event(1), event(2)], photons=[photon(1, 11, 21.0), photon(2, 12, 21.0)],
                truth_photons=[signal_truth_photon(0, 1, 1, 9, 1, 22.0),
                               signal_truth_photon(0, 1, 2, 9, 2, 22.0, hepmc_association_valid=0)],   # class kept, association missing
                links=[match_link(0, 1, 1, 5, 11, 9, 1), match_link(0, 1, 2, 5, 12, 9, 2)],
            )
            response = build_photon_response([path], system="pp", selection=SELECTION, truth_signal=truth_signal,
                                             event_weights={(0, 1, 1): 3.0})
            self.assertEqual(response.matrix[3, 2], 3.0)
            self.assertEqual(response.provenance["counts"]["dropped_no_weight"], 1)
            response = build_photon_response([path], system="pp", selection=SELECTION, truth_signal=truth_signal)
            self.assertEqual(response.fakes[2], 1.0)     # event 2: reco leader, no admitted truth photon


class HistogramWeightsAndLeakageTest(unittest.TestCase):
    def config(self):
        config = load_config(ROOT / "config" / "nominal.yaml")
        for system in ("pp", "auau"):
            config["systems"][system]["event"]["abs_vertex_z_max_cm"] = 1.0e9
            config["systems"][system]["event"].pop("centrality_percent", None)
        return config

    def test_leakage_keeps_one_candidate_per_truth_photon(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = write_tree_file(
                Path(temporary) / "sim.root", events=[event(1)],
                photons=[photon(1, 11, 21.0, ordinal=0), photon(1, 12, 23.0, ordinal=1)],
                truth_photons=[signal_truth_photon(0, 1, 1, 9, 1, 22.0)],
                links=[match_link(0, 1, 1, 5, 11, 9, 1, metric=0.01), match_link(0, 1, 1, 5, 12, 9, 1, metric=0.05)],
            )
            payload = make_histograms.build([(path, None)], system="pp", config=self.config(), thresholds_from="tree", truth_signal_only=True)
            counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
            self.assertEqual(counts[0].sum(), 1.0)                 # region A: one event-leading photon
            self.assertEqual(payload["abcd_event_leading_events_ptgamma"][0][1], 1)   # the 21 GeV best match, not the 23 GeV one

    def test_sample_weights_reach_the_histograms(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = write_tree_file(
                Path(temporary) / "auau_jet20.root",
                events=[event(1, centrality=12.0), event(2, centrality=12.0)],
                photons=[photon(1, 11, 21.0), photon(2, 12, 21.0)],
                truth_jets=[{"source_file_index": 0, "event_id_hi": 1, "event_id_lo": 1, "truth_jet_id_hi": 3, "truth_jet_id_lo": 1,
                             "truth_jet_radius": 0.4, "truth_jet_pt": 25.0, "truth_jet_eta": 0.0, "truth_jet_phi": 0.0},
                            {"source_file_index": 0, "event_id_hi": 1, "event_id_lo": 2, "truth_jet_id_hi": 3, "truth_jet_id_lo": 2,
                             "truth_jet_radius": 0.4, "truth_jet_pt": 40.0, "truth_jet_eta": 0.0, "truth_jet_phi": 0.0}],
            )
            payload = make_histograms.build([(path, "auau_jet20")], system="auau", config=self.config(), thresholds_from="tree",
                                            truth_signal_only=False, manifest=MANIFEST)
            expected = MANIFEST.samples["auau_jet20"].factor_value * MANIFEST.centrality.values["inclusivejet"][2]
            counts = np.asarray(payload["abcd_event_leading_counts_ptgamma"])
            self.assertAlmostEqual(counts[0].sum(), expected, places=15)          # event 2 (40 GeV) is outside the window
            self.assertAlmostEqual(np.asarray(payload["abcd_event_leading_sumw2_ptgamma"])[0].sum(), expected ** 2, places=20)
            self.assertEqual(payload["weights"]["per_input"][0]["dropped"]["outside_ownership_window"], 1)


if __name__ == "__main__":
    unittest.main()
