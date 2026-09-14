"""Bounded interface tests. Fixtures certify mechanics, not producer physics."""
from pathlib import Path
import hashlib
import math
import os
import shutil
import tempfile
import unittest
from unittest import mock

import awkward as ak
import numpy as np
import uproot

import build_photonjet_collaboration_tree as builder
import validate_photonjet_collaboration_tree as validator


def add_mbd(events):
    n = len(events["event_id_hi"])
    events["mbd_pmt_available"] = np.ones(n, dtype="int32")
    for key, values in {
        "mbd_pmt_id": list(range(128)), "mbd_pmt_arm": [0] * 64 + [1] * 64,
        "mbd_pmt_charge": [-1., 0.] + [2.] * 126,
        "mbd_pmt_charge_valid": [1] * 128,
    }.items():
        events[key] = ak.Array([values for _ in range(n)])


def event_fixture(complete=True):
    z = np.asarray([0., 29.999, 30., -30., 59.999, 60., -60., np.nan])
    n = len(z)
    rows = {name: np.zeros(n, dtype=dtype) for name, dtype in builder.LEGACY_EVENT_SCHEMA.items()
            if name not in ("source_file_index", "source_entry")}
    rows.update(event_id_hi=np.ones(n, dtype="uint64"),
                event_id_lo=np.arange(1, n + 1, dtype="uint64"),
                vertex_z=z, centrality=np.full(n, 25.), event_weight=np.full(n, 210.),
                source_occurrence_id_hi=np.full(n, 99, dtype="uint64"),
                source_occurrence_id_lo=np.full(n, 101, dtype="uint64"),
                truth_vertex_z=np.arange(n, dtype="float64") * 10,
                truth_mb_vertex_z=np.full(n, np.nan),
                mbd_total_charge=np.asarray([0., 3., np.nan, 10., 11., 12., 13., 14.]),
                nodes_ready=np.asarray([1] * 7 + [0], dtype="int32"))
    if complete:
        rows.update(reco_vertex_valid=np.asarray([1] * 7 + [0], dtype="int32"),
                    sample_weight_valid=np.ones(n, dtype="int32"),
                    vertex_weight_valid=np.asarray([1] * 7 + [0], dtype="int32"),
                    truth_denominator_complete=np.ones(n, dtype="int32"))
    return rows


def weight_fixture(events, system="pp"):
    n = len(events["event_id_hi"])
    rows = {name: np.ones(n) for name in builder.WEIGHT_FIELDS}
    rows.update(target_id_hi=events["event_id_hi"].copy(),
                target_id_lo=events["event_id_lo"].copy(),
                component_type=np.asarray(["event"] * n, dtype=object),
                application_count=np.ones(n, dtype="int32"))
    if system == "pp":
        rows.update(slice_weight=np.full(n, 3.), cross_section_weight=np.full(n, 3.),
                    vertex_weight=np.full(n, 2.), si_di_weight=np.full(n, 5.),
                    period_weight=np.full(n, 7.), exposure_weight=np.full(n, 7.),
                    final_weight=np.full(n, 210.))
    else:
        rows.update(vertex_weight=np.full(n, 2.), exposure_weight=np.full(n, 7.),
                    final_weight=np.full(n, 14.))
        events["event_weight"][:] = 14.
    # An undefined reco correction never kills the independent sample factor.
    rows["vertex_weight"][-1] = np.nan
    rows["final_weight"][-1] = np.nan
    events["event_weight"][-1] = np.nan
    return rows


def write_source(path, system="pp", simulation=True, complete=True, truth_cones=None, dominant=None,
                 source_occurrence=101, candidate_et=20., centrality=25., features_valid=True,
                 truth_vertices=None, jet_view_provenance=True,
                 truth_only_jet_event=False, schema13_witnesses=False):
    events = event_fixture(complete)
    events["source_occurrence_id_lo"][:] = source_occurrence
    events["centrality"][:] = centrality
    if truth_vertices is not None:
        events.update(truth_vertices)
    weights = weight_fixture(events, system)
    if complete:
        # A complete fixture certifies finite truth-vertex factors even for an
        # event without a reconstructed vertex. Raw/unknown-factor fixtures
        # above intentionally retain their NaN and invalid witnesses.
        events["vertex_weight_valid"][:] = 1
        weights["vertex_weight"][-1] = 2.
        weights["final_weight"][-1] = 210. if system == "pp" else 14.
        events["event_weight"][-1] = weights["final_weight"][-1]
        add_mbd(events)
    n = len(events["event_id_hi"])
    if schema13_witnesses:
        for name, dtype in builder.PPG12_EVENT_WITNESS_SCHEMA.items():
            value = (0.8 if name.endswith("_eta") else 42.5
                     if dtype.startswith("float") else 7)
            events[name] = np.full(n, value, dtype=dtype)

    def table(root, name, values):
        schema = {key: ("string" if array.dtype == object else str(array.dtype))
                  for key, array in values.items() if not isinstance(array, ak.Array)}
        schema.update({key: {**builder.MBD_ARRAY_SCHEMA, **builder.TRUTH_VERTEX_ARRAY_SCHEMA, **builder.TRUTH_JET_AVAILABILITY_SCHEMA}.get(key, "var * float64") for key, array in values.items()
                       if isinstance(array, ak.Array)})
        tree = root.mktree(builder.PREFIX + name, schema)
        if len(next(iter(values.values()))):
            tree.extend(values)

    def objects(branches, count, identity):
        ints = {"encounter_ordinal", "deterministic_order", "truth_signal_match_state",
                "truth_signal_match_barcode", "pass_state", "recoil_state", "photon_rank",
                "jet_rank", "wrong_photon_class", "wrong_recoil_class", "prompt_class",
                "source_role", "generator_barcode", "reco_type", "truth_type", "link_class",
                "truth_isolation_valid", "g4_photon_valid", "hepmc_association_valid",
                "analysis_signal_r03", "native_track_id", "native_vertex_id", "embedding_id"}
        ints.update(key for key in builder.CANDIDATE_INPUT_BRANCHES
                    if key.startswith("dominant_truth_") and not key.endswith("energy_contribution"))
        strings = {"input_identity", "subtraction_identity"}
        rows = {key: np.zeros(count, dtype=(object if key in strings else
                "uint64" if key.endswith(("_hi", "_lo")) or key == "quality_bitmask"
                else "int32" if key in ints else "float64"))
                for key in branches}
        if "event_id_hi" in rows:
            rows["event_id_hi"][:] = 1
            rows["event_id_lo"][:] = np.arange(1, count + 1)
        for key in rows:
            if key.endswith("_hi") and key != "event_id_hi":
                rows[key][:] = identity
            if key.endswith("_lo") and key != "event_id_lo":
                rows[key][:] = np.arange(1, count + 1)
        return rows

    with uproot.recreate(path) as root:
        table(root, "RJEventV1", events)
        table(root, "RJWeightComponentV1", weights)
        photons = objects(builder.CANDIDATE_INPUT_BRANCHES, 1, 2)
        photons["cluster_et"][:] = candidate_et
        if schema13_witnesses:
            photons["native_cluster_key"][:] = 12345
            photons["ppg12_tower_mask_state"][:] = 1
            photons["ppg12_source_eligible_state"][:] = 1
            photons["reference_preselection_state"][:] = 0
            photons["active_preselection_state"][:] = 1
        if dominant is not None:
            for key, value in dominant.items():
                photons[key][:] = value
        table(root, "RJPhotonCandidateV1", photons)
        dim = 11 if system == "pp" else 14
        features = ak.Array([[candidate_et, .1 if features_valid else np.nan] + [0.1] * (dim - 2)])
        table(root, "RJModelEvaluationV1", dict(candidate_id_hi=np.asarray([2], dtype="uint64"),
              candidate_id_lo=np.asarray([1], dtype="uint64"), raw_score=np.asarray([.9]),
              shower_definition_id=np.asarray(["H70"], dtype=object), ordered_input_witnesses=features))
        table(root, "RJShowerFeatureViewV1", dict(candidate_id_hi=np.asarray([2], dtype="uint64"),
              candidate_id_lo=np.asarray([1], dtype="uint64"), definition_name=np.asarray(["H70"], dtype=object),
              ordered_features=features, finite_feature_state=np.full(1, int(features_valid), dtype="int32")))
        iso = objects(builder.ISOLATION_INPUT_BRANCHES, 2, 2)
        iso["candidate_id_lo"][:] = 1
        iso["radius"][:] = [.3, .4]
        iso["threshold"][:] = 1
        iso["sideband_threshold"][:] = 3
        iso["pass_state"][:] = 1
        table(root, "RJIsolationWitnessV1", iso)
        jets = objects(builder.JET_INPUT_BRANCHES, 1, 3)
        jets["radius"][:] = .4
        jets["raw_pt"][:] = jets["corrected_pt"][:] = 10
        jets["phi"][:] = math.pi
        jets["input_identity"][:] = "towerinfo" if system == "pp" else "towerinfo_sub1"
        jets["subtraction_identity"][:] = "none" if system == "pp" else "SUB1"
        if not jet_view_provenance:
            del jets["input_identity"]
            del jets["subtraction_identity"]
        table(root, "RJJetV1", jets)
        pairs = objects(builder.PAIR_INPUT_BRANCHES, 1, 4)
        pairs["candidate_id_hi"][:] = 2
        pairs["jet_id_hi"][:] = 3
        pairs["delta_phi"][:] = math.pi
        pairs["xjgamma"][:] = 10 / candidate_et
        table(root, "RJPhotonJetPairV1", pairs)
        truth = objects(builder.TRUTH_PHOTON_INPUT_BRANCHES, n if simulation else 0, 5)
        truth["pt"][:] = 20
        truth["prompt_class"][:] = 1
        if schema13_witnesses:
            truth["source_role"][:] = 0
            truth["sample_source_role"][:] = 2
            truth["generator_occurrence_embedding_id"][:] = 4
        if truth_cones is not None:
            truth["truth_isolation_r03"][:] = truth_cones[0]
            truth["truth_isolation_r04"][:] = truth_cones[1]
            truth["truth_isolation_witness"][:] = truth_cones[0]
            truth["truth_isolation_valid"][:] = 1
            truth["g4_photon_valid"][:] = 1
            truth["hepmc_association_valid"][:] = 1
            truth["analysis_signal_r03"][:] = truth_cones[0] < 4
            truth["native_track_id"][:] = np.arange(n) + 10
            truth["embedding_id"][:] = 1
        table(root, "RJTruthPhotonV1", truth)
        truth_jet_count = (2 if truth_only_jet_event else 1) if simulation else 0
        tj = objects(builder.TRUTH_JET_INPUT_BRANCHES, truth_jet_count, 6)
        tj["pt"][:] = 10
        tj["radius"][:] = .4
        table(root, "RJTruthJetV1", tj)
        links = objects(builder.LINK_INPUT_BRANCHES, n if simulation else 0, 7)
        if simulation:
            links["reco_type"][0] = 1
            links["reco_id_hi"][:] = 0
            links["reco_id_lo"][:] = 0
            links["reco_id_hi"][0] = 2
            links["reco_id_lo"][0] = 1
            links["truth_type"][:] = 1
            links["truth_id_hi"][:] = 5
            links["link_class"][:] = 2
            links["link_class"][0] = 0
        table(root, "RJRecoTruthLinkV1", links)


class StorageContractTests(unittest.TestCase):
    def test_truth_inventory_empty_missing_and_contradictory_states(self):
        events = {"event_id_hi": np.zeros(2, dtype="uint64"),
                  "truth_denominator_complete": np.array([1, 0])}
        missing = builder._truth_inventory_columns(events)
        self.assertEqual(missing["truth_photon_capture_state"].tolist(), [-2, -2])
        for name in builder.TRUTH_INVENTORY_SCHEMA:
            events[name] = np.zeros(2, dtype="int32")
        events["truth_photon_capture_state"][:] = [1, 0]
        retained = builder._truth_inventory_columns(events)
        self.assertEqual(retained["truth_photon_capture_state"].tolist(), [1, 0])
        events["truth_photon_native_embedded_primary_pid22_count"][0] = 1
        with self.assertRaisesRegex(ValueError, "contradicts inventory"):
            builder._truth_inventory_columns(events)
        events["truth_denominator_complete"][0] = 0
        builder._truth_inventory_columns(events)
        events.pop("truth_photon_duplicate_track_count")
        with self.assertRaisesRegex(ValueError, "partial truth-photon"):
            builder._truth_inventory_columns(events)

    def test_schema13_ppg12_witnesses_survive_normalized_export(self):
        for system in ("pp", "auau"):
            with self.subTest(system=system), tempfile.TemporaryDirectory() as temp:
                source, output = (Path(temp) / name for name in ("source.root", "package.root"))
                write_source(source, system=system, schema13_witnesses=True)
                builder.build([source], output, system, require_complete_interface=True,
                              layout="normalized_v1")
                validator.validate([source], output, system, require_complete_interface=True)
                with uproot.open(source) as original, uproot.open(output) as package:
                    for name in builder.PPG12_EVENT_WITNESS_SCHEMA:
                        np.testing.assert_array_equal(
                            original[builder.PREFIX + "RJEventV1"][name].array(library="np"),
                            package["events"][name].array(library="np"))
                    for name, expected in (("native_cluster_key", 12345),
                                           ("ppg12_tower_mask_state", 1),
                                           ("ppg12_source_eligible_state", 1),
                                           ("reference_preselection_state", 0),
                                           ("active_preselection_state", 1)):
                        self.assertEqual(package["photons"][name].array(library="np").tolist(), [expected])
                    for name, expected in (("source_role", 0), ("sample_source_role", 2),
                                           ("generator_occurrence_embedding_id", 4)):
                        self.assertTrue(np.all(package["truthPhotons"][name].array(library="np") == expected))

    def test_per_view_response_preserves_candidates_and_labels_truth_misses(self):
        jets = {
            "event_id_hi": np.asarray([2, 2], dtype="uint64"),
            "event_id_lo": np.asarray([1, 1], dtype="uint64"),
            "jet_id_hi": np.asarray([8, 8], dtype="uint64"),
            "jet_id_lo": np.asarray([1, 2], dtype="uint64"),
            "radius": np.asarray([.4, .4]), "raw_pt": np.asarray([18., 30.]),
            "corrected_pt": np.asarray([18., 30.]), "eta": np.asarray([.06, .01]),
            "phi": np.asarray([3., 3.]), "mass": np.zeros(2), "area": np.ones(2),
            "quality_bitmask": np.zeros(2, dtype="uint64"),
            "deterministic_order": np.asarray([0, 0], dtype="int32"),
            "input_identity": np.asarray(
                ["towerinfo_sub1", "towerinfo_retower_unsubtracted"], dtype=object
            ),
            "subtraction_identity": np.asarray(["SUB1", "none"], dtype=object),
        }
        truth = {
            "event_id_hi": np.asarray([2, 2, 2], dtype="uint64"),
            "event_id_lo": np.asarray([1, 1, 2], dtype="uint64"),
            "truth_jet_id_hi": np.asarray([13, 13, 13], dtype="uint64"),
            "truth_jet_id_lo": np.asarray([1, 2, 3], dtype="uint64"),
            "radius": np.asarray([.4, .4, .4]), "pt": np.asarray([20., 15., 12.]),
            "eta": np.asarray([0., 1., 0.]), "phi": np.asarray([3., 3., 3.]),
        }
        links = {
            "link_id_hi": np.full(6, 14, dtype="uint64"),
            "link_id_lo": np.arange(1, 7, dtype="uint64"),
            "reco_type": np.asarray([2, 2, 2, 2, 0, 0], dtype="int32"),
            "reco_id_hi": np.asarray([8, 8, 8, 8, 0, 0], dtype="uint64"),
            "reco_id_lo": np.asarray([1, 2, 2, 1, 0, 0], dtype="uint64"),
            "truth_type": np.asarray([2, 2, 2, 0, 2, 2], dtype="int32"),
            "truth_id_hi": np.asarray([13, 13, 13, 0, 13, 13], dtype="uint64"),
            "truth_id_lo": np.asarray([1, 1, 1, 0, 2, 3], dtype="uint64"),
            "match_metric": np.asarray([.06, .01, .01, .06, .94, np.nan]),
            "link_class": np.asarray([5, 5, 0, 1, 2, 2], dtype="int32"),
        }
        reco_lookup = builder._local_identity_lookup(
            jets, "jet_id_hi", "jet_id_lo", "reco jet identity"
        )
        truth_lookup = builder._local_identity_lookup(
            truth, "truth_jet_id_hi", "truth_jet_id_lo", "truth jet identity"
        )
        rows = builder._jet_response_rows(
            0, jets, truth, links, reco_lookup, truth_lookup,
            {(2, 1): [("towerinfo_sub1", "SUB1", .4),
                      ("towerinfo_retower_unsubtracted", "none", .4)],
             (2, 2): [("towerinfo_sub1", "SUB1", .4)]},
        )
        candidates = [row for row in rows if row["jet_response_state"] == 2]
        finals = [row for row in rows if row["jet_response_state"] == 1]
        unsupported = [row for row in rows if row["jet_response_state"] == 0]
        self.assertEqual([(int(r["link_id_hi"]), int(r["link_id_lo"])) for r in candidates],
                         [(14, 1), (14, 2)])
        self.assertEqual([r["link_class"] for r in finals].count(0), 2)
        self.assertEqual([r["link_class"] for r in finals].count(1), 0)
        self.assertEqual([r["link_class"] for r in finals].count(2), 3)
        misses = [row for row in finals if row["link_class"] == 2]
        self.assertEqual({(r["jet_input_identity"], r["jet_subtraction_identity"])
                          for r in misses},
                         {("towerinfo_sub1", "SUB1"),
                          ("towerinfo_retower_unsubtracted", "none")})
        self.assertTrue(all(r["reco_type"] == 0 and r["reco_index"] == -1 for r in misses))
        self.assertEqual(len(unsupported), 0)
        event2_miss = [r for r in misses if int(r["event_id_lo"]) == 2]
        self.assertEqual(len(event2_miss), 1)
        self.assertEqual((event2_miss[0]["jet_input_identity"],
                          event2_miss[0]["jet_subtraction_identity"]),
                         ("towerinfo_sub1", "SUB1"))
        identities = [(int(r["link_id_hi"]), int(r["link_id_lo"])) for r in rows]
        self.assertEqual(len(identities), len(set(identities)))
        absent = builder._jet_response_rows(
            0, jets, truth, links, reco_lookup, truth_lookup,
            {(2, 1): [("towerinfo_sub1", "SUB1", .4),
                      ("towerinfo_retower_unsubtracted", "none", .4)]},
        )
        unsupported = [row for row in absent if row["jet_response_state"] == 0]
        self.assertEqual([(int(r["event_id_lo"]), r["jet_input_identity"],
                           r["jet_subtraction_identity"]) for r in unsupported],
                         [(2, "", "")])

    def test_jet_view_provenance_is_fail_closed_or_explicitly_unsupported(self):
        base = {"jet_id_hi": np.asarray([8], dtype="uint64")}
        inputs, subtractions, state = builder._jet_view_columns(base)
        self.assertEqual((inputs.tolist(), subtractions.tolist(), state.tolist()),
                         ([""], [""], [0]))
        with self.assertRaisesRegex(ValueError, "partial reconstruction-view provenance"):
            builder._jet_view_columns({**base, "input_identity": np.asarray(["towerinfo"])})
        with self.assertRaisesRegex(ValueError, "invalid reconstruction-view provenance"):
            builder._jet_view_columns({**base, "input_identity": np.asarray([""]),
                                       "subtraction_identity": np.asarray(["none"])})
        availability = {
            "event_id_hi": np.asarray([2, 2], dtype="uint64"),
            "event_id_lo": np.asarray([1, 2], dtype="uint64"),
            "reco_jet_view_radius_code": ak.Array([[4, 4], [4]]),
            "reco_jet_view_input_identity": ak.Array(
                [["towerinfo_sub1", "towerinfo_retower_unsubtracted"],
                 ["towerinfo_sub1"]]
            ),
            "reco_jet_view_subtraction_identity": ak.Array([["SUB1", "none"], ["SUB1"]]),
            "reco_jet_view_available": ak.Array([[1, 1], [1]]),
        }
        views, rows = builder._reco_jet_view_availability(availability)
        self.assertEqual(views[(2, 2)], [("towerinfo_sub1", "SUB1", .4)])
        self.assertEqual(len(rows), 3)
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "legacy.root", Path(temp) / "portable.root"
            write_source(source, jet_view_provenance=False)
            builder.build([source], output, "pp", layout="normalized_v1")
            report = validator.validate([source], output, "pp")
            self.assertIn("per_view_jet_response=UNSUPPORTED_MISSING_VIEW_AVAILABILITY_OR_LEGACY_PROVENANCE",
                          report)
            with uproot.open(output) as root:
                self.assertEqual(root["jets"]["jet_view_provenance_state"].array(
                    library="np").tolist(), [0])

    def test_truth_only_event_without_view_availability_is_not_vacuous_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "compact.root"
            # Event 1 proves this view/radius exists elsewhere in the run; it
            # does not prove that the empty reconstructed collection ran in event 2.
            write_source(source, truth_only_jet_event=True)
            builder.build([source], output, "pp", layout="normalized_v1")
            report = validator.validate([source], output, "pp")
            self.assertIn(
                "per_view_jet_response=UNSUPPORTED_MISSING_VIEW_AVAILABILITY_OR_LEGACY_PROVENANCE",
                report,
            )
            with uproot.open(output) as root:
                rows = root["jetResponseLinks"].arrays(library="np")
                unsupported = rows["jet_response_state"] == 0
                self.assertEqual(int(np.count_nonzero(unsupported)), 2)
                self.assertEqual(rows["event_id_lo"][unsupported].tolist(), [1, 2])
                self.assertEqual(rows["link_class"][unsupported].tolist(), [2, 2])
                self.assertEqual(rows["jet_input_identity"][unsupported].tolist(), ["", ""])
                self.assertEqual(rows["jet_subtraction_identity"][unsupported].tolist(), ["", ""])

    def test_input_aliases_and_output_collision_fail_before_write(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            source, output = folder / "source.root", folder / "output.root"
            write_source(source)
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            symbolic, hard = folder / "symbolic.root", folder / "hard.root"
            symbolic.symlink_to(source)
            os.link(source, hard)
            for duplicate in (source, symbolic, hard):
                with self.subTest(duplicate=duplicate):
                    with self.assertRaisesRegex(ValueError, "duplicate input"):
                        builder.build([source, duplicate], output, "pp")
                    self.assertFalse(output.exists())
                    with self.assertRaisesRegex(ValueError, "duplicate input"):
                        validator.validate([source, duplicate], output, "pp")
                    with self.assertRaisesRegex(ValueError, "output aliases"):
                        builder.build([source], duplicate, "pp")
                    self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), digest)
            copied = folder / "copied.root"
            shutil.copyfile(source, copied)
            with self.assertRaisesRegex(ValueError, "duplicate ROOT file UUID"):
                builder.build([source, copied], output, "pp")
            self.assertFalse(output.exists())

    def test_normalized_roundtrip_matches_expanded_payloads(self):
        for system in ("pp", "auau"):
            with self.subTest(system=system), tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                source, expanded, compact = (folder / name for name in
                                             ("source.root", "expanded.root", "compact.root"))
                write_source(source, system=system, truth_cones=(1.5, 4.5))
                old = builder.build([source], expanded, system, require_complete_interface=True)
                new = builder.build([source], compact, system, require_complete_interface=True,
                                    layout="normalized_v1")
                self.assertEqual(new, {k: v for k, v in old.items() if k != "eventTree"})
                report = validator.validate([source], compact, system, require_complete_interface=True)
                self.assertIn("normalized_single_payload_layout=PASS", report)
                self.assertIn("normalized_pair_joins_and_kinematics=PASS", report)
                with uproot.open(expanded) as before, uproot.open(compact) as after:
                    self.assertNotIn("eventTree", after)
                    for tree in new:
                        for name in builder.output_schemas("normalized_v1")[tree]:
                            old_tree = ("eventTree" if tree == "events" and
                                        name not in builder.EVENT_SCHEMA else tree)
                            a = after[tree][name].array(library="ak")
                            b = before[old_tree][name].array(library="ak")
                            self.assertEqual(ak.to_list(ak.num(a, axis=0)), ak.to_list(ak.num(b, axis=0)))
                            np.testing.assert_array_equal(ak.to_numpy(ak.flatten(a, axis=None)),
                                                          ak.to_numpy(ak.flatten(b, axis=None)))
                    self.assertEqual(np.count_nonzero(after["recoTruthLinks"]["reco_type"].array(library="np") == 0), 7)
                self.assertLess(compact.stat().st_size, expanded.stat().st_size)
                relocated = folder / "relocated_source.root"
                shutil.copyfile(source, relocated)
                validator.validate([relocated], compact, system, require_complete_interface=True)

    def test_normalized_unknown_context_and_guard_features_survive(self):
        cases = ((20., np.nan, True), (20., 85., True), (10., 25., False), (20., 25., False))
        for et, cent, valid in cases:
            for supplied_model in (False, True):
                with self.subTest(et=et, centrality=cent, valid=valid, model=supplied_model), tempfile.TemporaryDirectory() as temp:
                    source, output = Path(temp) / "source.root", Path(temp) / "compact.root"
                    write_source(source, system="auau", candidate_et=et, centrality=cent, features_valid=valid)
                    kwargs = dict(model_path=Path("fixture-model.root"), model_input_count=14) if supplied_model else {}
                    evaluator = (lambda values: .8) if supplied_model else None
                    with mock.patch.object(builder, "_make_rbdt_evaluator", return_value=evaluator), \
                         mock.patch.object(validator, "_make_evaluator", return_value=evaluator):
                        counts = builder.build([source], output, "auau", layout="normalized_v1", **kwargs)
                        report = validator.validate([source], output, "auau", **kwargs)
                    self.assertEqual(counts["photons"], 1)
                    self.assertEqual(counts["photonJets"], 1)
                    self.assertIn("scored_candidates=0; unscored_candidates=1", report)
                    with uproot.open(output) as root:
                        rows = root["photons"].arrays(library="np")
                        self.assertTrue(np.isnan(rows["bdt_score"][0]))
                        self.assertEqual(rows["bdt_is_not_tight"][0], -1)
                        self.assertEqual(rows["bdt_input_count"][0], 14)
                        self.assertEqual(rows["photon_et"][0], et)

    def test_normalized_pair_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "compact.root"
            write_source(source)
            builder.build([source], output, "pp", layout="normalized_v1")
            with uproot.open(output) as root:
                original = root["photonJets"].arrays(library="np")
            for field, value in (("candidate_id_lo", 999), ("photon_index", 3),
                                  ("source_entry", 7), ("xjgamma", 1.5)):
                with self.subTest(field=field):
                    rows = {k: v.copy() for k, v in original.items()}
                    rows[field][0] = value
                    with uproot.update(output) as root:
                        root.mktree("photonJets", builder.output_schemas("normalized_v1")["photonJets"]).extend(rows)
                    with self.assertRaisesRegex(ValueError, "normalized pair"):
                        validator.validate([source], output, "pp")

    def test_normalized_supplied_model_is_used_for_supported_rows(self):
        for system, dimension in (("pp", 11), ("auau", 14)):
            with self.subTest(system=system), tempfile.TemporaryDirectory() as temp:
                source, output = Path(temp) / "source.root", Path(temp) / "compact.root"
                write_source(source, system=system)
                kwargs = dict(model_path=Path("fixture-model.root"), model_input_count=dimension)
                with mock.patch.object(builder, "_make_rbdt_evaluator", return_value=lambda v: .81), \
                     mock.patch.object(validator, "_make_evaluator", return_value=lambda v: .81):
                    builder.build([source], output, system, layout="normalized_v1", **kwargs)
                    report = validator.validate([source], output, system, **kwargs)
                self.assertIn("scored_candidates=1; unscored_candidates=0", report)
                with uproot.open(output) as root:
                    self.assertEqual(root["photons"]["bdt_evaluation_state"].array(library="np").tolist(), [1])
                with uproot.update(output) as root:
                    root["source_files"] = "[]"
                with self.assertRaisesRegex(ValueError, "source-file locator"):
                    validator.validate([source], output, system)

    def test_legacy_filtered_guard_matches_dominant_readback(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "expanded.root"
            write_source(source, candidate_et=10., features_valid=False)
            with mock.patch.object(builder, "_make_rbdt_evaluator", return_value=lambda v: .8), \
                 mock.patch.object(validator, "_make_evaluator", return_value=lambda v: .8):
                kwargs = dict(model_path=Path("fixture-model.root"), model_input_count=11)
                counts = builder.build([source], output, "pp", **kwargs)
                report = validator.validate([source], output, "pp", **kwargs)
            self.assertEqual(counts["photons"], 0)
            self.assertIn("dominant_primary_source_readback=PASS", report)


class InterfaceTests(unittest.TestCase):
    def test_native_vertex_arrays_survive_both_package_layouts(self):
        vertices = {
            "truth_vertex_capture_state": np.ones(8, dtype="int32"),
            "truth_primary_vertex_id": np.full(8, 2, dtype="int32"),
            "truth_vertex_id": ak.Array([[1,2]]*8),
            "truth_vertex_embedding_id": ak.Array([[0,1]]*8),
            "truth_vertex_embedding_valid": ak.Array([[1,1]]*8),
            "truth_vertex_z_valid": ak.Array([[1,1]]*8),
            "truth_vertex_z_values": ak.Array([[-12.,18.],[-12.,52.]]*4),
            "truth_jet_capture_state": np.ones(8, dtype="int32"),
            "truth_jet_radius_code": ak.Array([[2,3,4]]*8),
            "truth_jet_container_valid": ak.Array([[1,0,1],[0,1,0]]*4),
        }
        with tempfile.TemporaryDirectory() as directory:
            for system in ("pp", "auau"):
                source = Path(directory)/f"{system}.root"
                write_source(source, system, truth_vertices=vertices)
                for layout in builder.LAYOUTS:
                    output = Path(directory)/f"{system}-{layout}.root"
                    builder.build([source], output, system, layout=layout)
                    validator.validate([source], output, system, None)
                    with uproot.open(output) as root:
                        tree = root["events" if layout=="normalized_v1" else "eventTree"]
                        names = list(builder.TRUTH_VERTEX_ARRAY_SCHEMA) + list(builder.TRUTH_JET_AVAILABILITY_SCHEMA)
                        actual = tree.arrays(names, library="ak")
                        for name in names:
                            self.assertEqual(ak.to_list(actual[name]), ak.to_list(vertices[name]))

    def test_truth_jet_availability_rejects_partial_and_contradictory_witnesses(self):
        with self.assertRaises(ValueError):
            builder._truth_jet_availability({"event_id_hi": np.ones(1, dtype="uint64"),
                "truth_jet_capture_state": np.ones(1, dtype="int32")})
        for states, keys, valid in (([0], [[4]], [[1]]), ([1], [[4,4]], [[1,1]]),
                                   ([1], [[4]], [[2]]), ([1], [[4]], [[]])):
            with self.assertRaises(ValueError):
                builder._truth_jet_availability({"event_id_hi": np.ones(1, dtype="uint64"),
                    "truth_jet_capture_state": np.array(states),
                    "truth_jet_radius_code": ak.Array(keys),
                    "truth_jet_container_valid": ak.Array(valid)})

    def setUp(self):
        self.events = event_fixture()
        self.weights = weight_fixture(self.events)

    def project(self, **kwargs):
        return builder._interface_columns(self.events, self.weights, "pp", **kwargs)

    def test_independent_factors_and_aliases(self):
        rows = self.project()
        np.testing.assert_array_equal(rows["sample_weight"], [3.] * 8)
        np.testing.assert_array_equal(rows["event_weight_without_vertex"], [105.] * 8)
        self.assertEqual(rows["weight_identity_valid"].tolist(), [1] * 7 + [0])
        self.assertTrue(np.isnan(rows["vertex_weight"][-1]))

    def test_strict_requires_vertex_certification_and_weight_identity(self):
        with self.assertRaisesRegex(ValueError, "vertex_weight_valid"):
            self.project(require_complete=True)
        self.events["vertex_weight_valid"][:] = 1
        self.weights["vertex_weight"][-1] = 2.
        self.weights["final_weight"][-1] = 210.
        self.events["event_weight"][-1] = 210.
        self.project(require_complete=True)
        for state in (0, -1):
            self.events["vertex_weight_valid"][0] = state
            with self.assertRaisesRegex(ValueError, "vertex_weight_valid"):
                self.project(require_complete=True)
        self.events["vertex_weight_valid"][0] = 1
        self.events["event_weight"][0] = np.nan
        with self.assertRaisesRegex(ValueError, "weight_identity_valid"):
            self.project(require_complete=True)

    def test_auau_existing_centrality_factor(self):
        w = weight_fixture(self.events, "auau")
        rows = builder._interface_columns(self.events, w, "auau")
        np.testing.assert_array_equal(rows["event_weight_without_vertex"], [7.] * 8)
        np.testing.assert_array_equal(rows["sample_weight"], [1.] * 8)

    def test_corrected_early_return_keeps_original_weight_witness(self):
        self.events["legacy_event_weight"] = np.ones(8)
        rows = self.project()
        np.testing.assert_array_equal(rows["legacy_event_weight"], np.ones(8))
        self.assertTrue(np.isnan(self.events["event_weight"][-1]))
        self.assertEqual(rows["sample_weight"][-1], 3.)
        self.assertEqual(rows["legacy_event_weight"][-1], 1.)
        del self.events["legacy_event_weight"]
        np.testing.assert_equal(self.project()["legacy_event_weight"], self.events["event_weight"])

    def test_strict_z_boundaries_and_missing_vertex(self):
        rows = self.project()
        self.assertEqual(rows["reco_vertex_abs_lt_30"].tolist(), [1, 1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(rows["reco_vertex_abs_lt_60"].tolist(), [1, 1, 1, 1, 1, 0, 0, 0])
        self.assertTrue(np.isnan(rows["reco_vertex_z"][-1]))

    def test_legacy_reset_values_do_not_certify_denominator(self):
        events = event_fixture(False)
        weights = weight_fixture(events)
        rows = builder._interface_columns(events, weights, "pp")
        self.assertTrue(np.all(rows["reco_vertex_valid"] == -1))
        self.assertTrue(np.all(rows["truth_denominator_complete"] == -1))
        with self.assertRaisesRegex(ValueError, "complete interface unavailable"):
            builder._interface_columns(events, weights, "pp", require_complete=True)

    def test_absent_weights_are_nan_not_unity(self):
        rows = builder._interface_columns(event_fixture(False), None, "pp")
        self.assertTrue(np.all(np.isnan(rows["sample_weight"])))
        self.assertTrue(np.all(rows["weight_components_recorded"] == 0))

    def test_duplicate_or_orphan_weights_fail(self):
        for value, message in ((1, "duplicate"), (90, "orphan")):
            with self.subTest(value=value):
                self.weights["target_id_lo"][-1] = value
                with self.assertRaisesRegex(ValueError, message):
                    self.project()

    def test_shuffled_join_and_non_event_component(self):
        indices = np.arange(7, -1, -1)
        self.weights = {key: values[indices] for key, values in self.weights.items()}
        rows = self.project()
        self.assertTrue(np.isnan(rows["vertex_weight"][-1]))
        self.weights = {key: np.concatenate([values, values[:1]]) for key, values in self.weights.items()}
        self.weights["component_type"][-1] = "photon"
        self.project()

    def test_missing_partial_weights_fail(self):
        del self.weights["exposure_weight"]
        with self.assertRaisesRegex(ValueError, "incomplete RJWeight"):
            self.project()

    def test_missing_join_row_fails(self):
        self.weights = {key: value[:-1] for key, value in self.weights.items()}
        with self.assertRaisesRegex(ValueError, "exactly one"):
            self.project()

    def test_numeric_mismatch_fails(self):
        self.events["event_weight"][0] = 105.
        with self.assertRaisesRegex(ValueError, "multiplication identity"):
            self.project()

    def test_alias_mismatch_fails(self):
        self.weights["cross_section_weight"][0] = 99.
        with self.assertRaisesRegex(ValueError, "alias mismatch"):
            self.project()

    def test_infinite_factor_is_not_valid(self):
        self.weights["vertex_weight"][0] = np.inf
        self.events["vertex_weight_valid"][0] = 0
        self.assertEqual(self.project()["weight_identity_valid"][0], 0)
        self.events["vertex_weight_valid"][0] = 1
        with self.assertRaisesRegex(ValueError, "requires finite"):
            self.project()

    def test_mbd_observed_zero_versus_absent(self):
        rows = self.project()
        self.assertEqual(rows["mbd_total_charge_finite"][:3].tolist(), [1, 1, 0])
        del self.events["mbd_total_charge"]
        self.assertTrue(np.all(self.project()["mbd_total_charge_finite"] == 0))

    def test_bad_witness_fails(self):
        self.events["reco_vertex_valid"][0] = 2
        with self.assertRaisesRegex(ValueError, "tri-state"):
            self.project()

    def test_per_pmt_charge_identity_order_and_unclipped_sums(self):
        add_mbd(self.events)
        scalars, arrays = builder._mbd_columns(self.events)
        self.assertEqual(ak.to_list(arrays["mbd_pmt_id"][0]), list(range(128)))
        self.assertEqual(scalars["mbd_pmt_south_charge"][0], 123.)
        self.assertEqual(scalars["mbd_pmt_north_charge"][0], 128.)
        self.assertEqual(scalars["mbd_pmt_sum_charge"][0], 251.)
        self.assertEqual(ak.to_list(arrays["mbd_pmt_charge"][0])[:2], [-1., 0.])

    def test_per_pmt_missing_invalid_and_duplicate(self):
        scalars, arrays = builder._mbd_columns(self.events)
        self.assertEqual(ak.to_list(ak.num(arrays["mbd_pmt_id"])), [0] * 8)
        self.assertTrue(np.all(scalars["mbd_pmt_available"] == -1))
        add_mbd(self.events)
        valid = ak.to_list(self.events["mbd_pmt_charge_valid"])
        valid[0][4] = 0
        self.events["mbd_pmt_charge_valid"] = ak.Array(valid)
        self.assertTrue(np.isnan(builder._mbd_columns(self.events)[0]["mbd_pmt_sum_charge"][0]))
        ids = ak.to_list(self.events["mbd_pmt_id"])
        ids[0][4] = 3
        self.events["mbd_pmt_id"] = ak.Array(ids)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            builder._mbd_columns(self.events)

    def test_per_pmt_missing_geometry_uses_documented_channel_map(self):
        add_mbd(self.events)
        original_charges = ak.to_list(self.events["mbd_pmt_charge"])
        self.events["mbd_pmt_arm"] = ak.Array([[-1] * 128] * 8)
        scalars, arrays = builder._mbd_columns(self.events)
        self.assertEqual(ak.to_list(arrays["mbd_pmt_arm"][0]), [0] * 64 + [1] * 64)
        self.assertEqual(scalars["mbd_pmt_arm_from_channel_map"].tolist(), [1] * 8)
        self.assertEqual(scalars["mbd_pmt_south_charge"][0], 123.)
        self.assertEqual(scalars["mbd_pmt_north_charge"][0], 128.)
        self.assertEqual(ak.to_list(arrays["mbd_pmt_charge"]), original_charges)
        self.assertEqual(ak.to_list(self.events["mbd_pmt_arm"][0]), [-1] * 128)
        # A conflicting observed arm is retained and prevents certified sums.
        arms = [[0] * 128] * 8
        self.events["mbd_pmt_arm"] = ak.Array(arms)
        scalars, arrays = builder._mbd_columns(self.events)
        self.assertEqual(scalars["mbd_pmt_arm_from_channel_map"].tolist(), [0] * 8)
        self.assertEqual(scalars["mbd_pmt_sums_valid"].tolist(), [0] * 8)

    def test_complete_eight_tree_roundtrip_and_truth_misses(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "part.root"
            write_source(source)
            second = Path(temp) / "second_source.root"
            write_source(second, source_occurrence=102)
            counts = builder.build([source, second], output, "pp",
                                   require_complete_interface=True)
            self.assertEqual(counts, dict(events=16, eventTree=16, photons=2, jets=2,
                photonJets=2, truthPhotons=16, truthJets=2, recoTruthLinks=16,
                jetResponseLinks=2, jetViewAvailability=0, triggerScalers=0, triggerRunInfo=0))
            report = validator.validate([source, second], output, "pp",
                                        require_complete_interface=True)
            self.assertIn("collaborator_interface_v2_source_and_event_joins=PASS", report)
            with uproot.open(output) as root:
                np.testing.assert_array_equal(root["eventTree"]["nphotons"].array(library="np"),
                                              [1, 0, 0, 0, 0, 0, 0, 0] * 2)
                self.assertEqual(np.count_nonzero(root["recoTruthLinks"]["reco_type"].array(library="np") == 0), 14)

    def test_dual_truth_cones_roundtrip_and_missed_photons(self):
        for system in ("pp", "auau"):
            with self.subTest(system=system), tempfile.TemporaryDirectory() as temp:
                source, output = Path(temp) / "source.root", Path(temp) / "part.root"
                write_source(source, system=system, truth_cones=(1.5, 4.5))
                builder.build([source], output, system, require_complete_interface=True)
                report = validator.validate([source], output, system, require_complete_interface=True)
                self.assertIn("truth_dual_cone_source_readback=PASS", report)
                with uproot.open(output) as root:
                    rows = root["truthPhotons"].arrays(library="np")
                    self.assertEqual(len(rows["truth_photon_pt"]), 8)
                    np.testing.assert_array_equal(rows["truth_isolation_r03"], [1.5] * 8)
                    np.testing.assert_array_equal(rows["truth_isolation_r04"], [4.5] * 8)
                    self.assertTrue(np.all(rows["truth_isolation_r03"] < 4))
                    self.assertFalse(np.any(rows["truth_isolation_r04"] < 4))
                    np.testing.assert_array_equal(rows["native_track_id"], np.arange(8) + 10)
                    self.assertEqual(np.count_nonzero(root["recoTruthLinks"]["reco_type"].array(library="np") == 0), 7)
                rows["truth_isolation_r04"][0] = rows["truth_isolation_r03"][0]
                with uproot.update(output) as root:
                    root.mktree("truthPhotons", builder.TRUTH_PHOTON_SCHEMA).extend(rows)
                with uproot.open(output) as root:
                    with self.assertRaisesRegex(ValueError, "truth_isolation_r04"):
                        validator.validate_truth_witnesses([source], root)

    def test_dominant_primary_roundtrip_and_unavailable_are_distinct(self):
        expected = dict(dominant_truth_state=2, dominant_truth_evaluator_mode=1,
                        dominant_truth_track_id=123, dominant_truth_pid=211,
                        dominant_truth_barcode=42, dominant_truth_embedding_id=3,
                        dominant_truth_energy_contribution=17.5)
        for system in ("pp", "auau"):
            with self.subTest(system=system), tempfile.TemporaryDirectory() as temp:
                source, output = Path(temp) / "source.root", Path(temp) / "part.root"
                write_source(source, system=system, dominant=expected)
                builder.build([source], output, system, require_complete_interface=True)
                with uproot.open(output) as root:
                    validator.validate_dominant_truth_witnesses([source], root, require_complete=True)
                    values = root["photons"].arrays(library="np")
                    for name, value in expected.items():
                        np.testing.assert_array_equal(values[name], [value])
                values["dominant_truth_track_id"][0] = 999
                with uproot.update(output) as root:
                    root.mktree("photons", {**builder.EVENT_SCHEMA, **builder.PHOTON_FIELDS}).extend(values)
                with uproot.open(output) as root:
                    with self.assertRaisesRegex(ValueError, "dominant_truth_track_id"):
                        validator.validate_dominant_truth_witnesses([source], root)
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "part.root"
            write_source(source)  # explicit evaluator-unavailable fixture
            builder.build([source], output, "pp", require_complete_interface=True)
            with uproot.open(output) as root:
                with self.assertRaisesRegex(ValueError, "not a background label"):
                    validator.validate_dominant_truth_witnesses([source], root, require_complete=True)

    def test_physical_event_validity_is_not_inferred_from_worker_counter(self):
        events = event_fixture()
        n = len(events["event_sequence"])
        events["physical_event_sequence"] = np.arange(n, dtype="int64") + 45678
        events["physical_event_sequence_valid"] = np.ones(n, dtype="int32")
        observed = builder._event_columns(events, 0)
        np.testing.assert_array_equal(observed["physical_event_sequence"], np.arange(n) + 45678)
        np.testing.assert_array_equal(observed["physical_event_sequence_valid"], np.ones(n))
        del events["physical_event_sequence_valid"]
        observed = builder._event_columns(events, 0)
        np.testing.assert_array_equal(observed["physical_event_sequence_valid"], np.full(n, -1))

    def test_mbd_required_rejects_old_source(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "part.root"
            write_source(source, complete=False)
            with self.assertRaisesRegex(ValueError, "per-PMT MBD unavailable"):
                builder.build([source], output, "pp", require_mbd_pmt=True)

    def test_truth_link_target_type_event_source_and_index(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "part.root"
            write_source(source)
            second = Path(temp) / "second_source.root"
            write_source(second, source_occurrence=102)
            builder.build([source, second], output, "pp", require_complete_interface=True)
            with uproot.open(output) as root:
                validator._validate_truth_link_targets(root)
                original = root["recoTruthLinks"].arrays(library="np")
            selected = int(np.flatnonzero(original["truth_type"] == 1)[0])
            for field, value, message in (
                ("truth_index", 999, "local index differs"),
                ("truth_type", 2, "typed truth target"),
                ("truth_type", 99, "unknown truth target type"),
                ("event_id_lo", 999, "target event or local index differs"),
                ("source_file_index", 999, "typed .* target"),
                ("truth_type", 0, "null local-index"),
                ("reco_index", 999, "local index differs"),
            ):
                with self.subTest(field=field, value=value):
                    changed = {name: values.copy() for name, values in original.items()}
                    changed[field][selected] = value
                    with uproot.update(output) as root:
                        root.mktree("recoTruthLinks", {name: values.dtype for name, values in changed.items()}).extend(changed)
                    with uproot.open(output) as root:
                        with self.assertRaisesRegex(ValueError, message):
                            validator._validate_truth_link_targets(root)

    def test_mbd_same_event_readback_detects_corruption(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source.root", Path(temp) / "part.root"
            write_source(source)
            builder.build([source], output, "pp")
            with uproot.open(output) as root:
                rows = root["eventTree"].arrays(library="ak")
            values = {name: rows[name] for name in
                      {**builder.EVENT_SCHEMA, **builder.EVENT_ARRAY_SCHEMA, **builder.MBD_ARRAY_SCHEMA,
                       **builder.centrality_replay.ARRAYS, **builder.TRUTH_VERTEX_ARRAY_SCHEMA, **builder.TRUTH_JET_AVAILABILITY_SCHEMA}}
            charge = ak.to_list(values["mbd_pmt_charge"])
            charge[0][3] += 1.
            values["mbd_pmt_charge"] = ak.Array(charge)
            with uproot.update(output) as root:
                root.mktree("eventTree", {**builder.EVENT_SCHEMA, **builder.EVENT_ARRAY_SCHEMA,
                                          **builder.MBD_ARRAY_SCHEMA, **builder.centrality_replay.ARRAYS,
                                          **builder.TRUTH_VERTEX_ARRAY_SCHEMA, **builder.TRUTH_JET_AVAILABILITY_SCHEMA}).extend(values)
            with uproot.open(output) as root:
                with self.assertRaisesRegex(ValueError, "MBD per-PMT source readback mismatch"):
                    validator.validate_interface([source], root, "pp")


if __name__ == "__main__":
    unittest.main()
