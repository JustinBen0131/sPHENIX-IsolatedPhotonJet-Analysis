from __future__ import annotations

from copy import deepcopy
from contextlib import redirect_stdout
import fcntl
import importlib.util
import io
import json
import shutil
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import numpy as np
import uproot

from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection, recoil_histogram, write_recoil_skim
from photonjet.analysis.response import Category
import photonjet.analysis.response_builder as response_builder_module
from photonjet.analysis.response_builder import (
    ResponseBuildConfig,
    build_response,
    write_response_artifacts,
)
from photonjet.cli import main as cli_main
from photonjet.plotting.contract import compile_plot_contract, write_histogram_receipt
from photonjet.provenance import write_json


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
GENERATOR_SPEC = importlib.util.spec_from_file_location(
    "photonjet_response_fixture_generator",
    ROOT / "tools" / "generate_engineering_fixtures.py",
)
if GENERATOR_SPEC is None or GENERATOR_SPEC.loader is None:
    raise RuntimeError("engineering-fixture generator cannot be loaded")
GENERATOR = importlib.util.module_from_spec(GENERATOR_SPEC)
GENERATOR_SPEC.loader.exec_module(GENERATOR)


class ResponseBuilderTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)

    @staticmethod
    def _specification(system: str) -> tuple[dict, dict]:
        seed = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(
                encoding="utf-8"
            )
        )
        contract = GENERATOR.load_contract()
        specification = deepcopy(
            GENERATOR.materialize_records(seed, contract)["systems"][system]
        )
        return specification, contract

    def _write(self, name: str, system: str, specification: dict, contract: dict) -> Path:
        output = self.output / name
        GENERATOR.write_fixture(
            output=output,
            system=system,
            specification=specification,
            contract=contract,
        )
        return output

    def test_both_systems_build_both_response_dimensions(self) -> None:
        for system in ("pp", "auau"):
            source = FIXTURES / f"photonjet_trees_{system}.root"
            for dimension in ("1D", "2D"):
                with self.subTest(system=system, dimension=dimension):
                    bundle = build_response(
                        [source],
                        ResponseBuildConfig(system=system, dimension=dimension),
                    )
                    self.assertEqual(bundle.system, system)
                    self.assertEqual(bundle.dimension, dimension)
                    self.assertTrue(
                        bundle.conservation()[
                            "closes_rtol_1e-10_atol_1e-10"
                        ]
                    )

    def test_auau_fixture_builds_matched_and_combinatoric_pairs(self) -> None:
        bundle = build_response(
            [FIXTURES / "photonjet_trees_auau.root"],
            ResponseBuildConfig(system="auau", dimension="2D"),
        )
        self.assertEqual(bundle.matrix.shape, (126, 105))
        self.assertEqual(float(bundle.matrix.sum()), 2.0)
        self.assertEqual(float(bundle.fakes.sum()), 2.0)
        self.assertEqual(float(bundle.misses.sum()), 0.0)
        self.assertEqual(float(bundle.reco.sum()), 4.0)
        self.assertEqual(float(bundle.truth.sum()), 2.0)
        self.assertEqual(
            float(bundle.fake_causes[Category.COMBINATORIC.value].sum()),
            2.0,
        )
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])
        self.assertEqual(
            bundle.status,
            "TREE_CONSTRUCTION_COMPLETE__GOLDEN_EQUIVALENCE_NOT_RUN",
        )
        self.assertEqual(
            bundle.unresolved["full_statistics_response_equivalence"],
            "NOT_RUN",
        )

    def test_pp_region_a_empty_reco_becomes_truth_misses(self) -> None:
        bundle = build_response(
            [FIXTURES / "photonjet_trees_pp.root"],
            ResponseBuildConfig(system="pp", dimension="1D"),
        )
        self.assertEqual(bundle.matrix.shape, (6, 5))
        self.assertEqual(float(bundle.reco.sum()), 0.0)
        self.assertEqual(float(bundle.truth.sum()), 2.0)
        self.assertEqual(float(bundle.misses.sum()), 2.0)
        self.assertEqual(
            bundle.provenance["observed"]["category_counts"],
            {Category.PHOTON_RECO_MISS.value: 2},
        )
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_event_weights_propagate_to_sumw2(self) -> None:
        specification, contract = self._specification("auau")
        for rows in specification["trees"].values():
            for row in rows:
                if row.get("event_id_lo") == 1001 and "event_weight" in row:
                    row["event_weight"] = 2.0
        source = self._write("weighted.root", "auau", specification, contract)
        bundle = build_response(
            [source], ResponseBuildConfig(system="auau", dimension="1D")
        )
        self.assertEqual(float(bundle.matrix.sum()), 3.0)
        self.assertEqual(float(bundle.matrix_sumw2.sum()), 5.0)
        self.assertEqual(float(bundle.fakes.sum()), 3.0)
        self.assertEqual(float(bundle.fakes_sumw2.sum()), 5.0)
        self.assertEqual(float(bundle.reco_sumw2.sum()), 10.0)
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_retained_terminal_event_is_excluded_from_physics_response(self) -> None:
        specification, contract = self._specification("auau")
        retained_event = 1001
        for rows in specification["trees"].values():
            for row in rows:
                if (
                    row.get("event_id_lo") == retained_event
                    and "terminal_status" in row
                ):
                    row["terminal_status"] = 2
        source = self._write("retained_terminal.root", "auau", specification, contract)
        bundle = build_response(
            [source], ResponseBuildConfig(system="auau", dimension="1D")
        )
        self.assertEqual(float(bundle.matrix.sum()), 1.0)
        self.assertEqual(float(bundle.fakes.sum()), 1.0)
        self.assertEqual(bundle.provenance["observed"]["input_events"], 2)
        self.assertEqual(bundle.provenance["observed"]["accepted_events"], 1)
        self.assertEqual(
            bundle.provenance["observed"]["retained_terminal_events_excluded"],
            1,
        )
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_terminal_event_exclusion_agrees_across_all_offline_products(self) -> None:
        """Keeping a rejected row must equal removing it before analysis.

        The oracle is a physical subset of the same input, not a second copy
        of the selection predicate. Diagnostic rows retain their payload.
        """
        for system in ("pp", "auau"):
            specification, contract = self._specification(system)
            accepted_only = deepcopy(specification)
            for name, rows in accepted_only["trees"].items():
                accepted_only["trees"][name] = [
                    row for row in rows if row.get("event_id_lo") != 1001
                ]
            for rows in specification["trees"].values():
                for row in rows:
                    if row.get("event_id_lo") == 1001 and "terminal_status" in row:
                        row["terminal_status"] = 2
            retained = self._write(f"{system}_retained.root", system, specification, contract)
            subset = self._write(f"{system}_accepted.root", system, accepted_only, contract)
            with uproot.open(retained) as root:
                self.assertEqual(root["events"].num_entries, 2)
                self.assertEqual(root["photons"].num_entries, 4)

            for definition in ("bounded", "complement"):
                for region in ("inclusive", "A", "B", "C", "D"):
                    with self.subTest(system=system, definition=definition, region=region):
                        selection = RecoilSelection(region=region, non_tight_definition=definition)
                        actual = recoil_histogram([retained], selection)
                        expected = recoil_histogram([subset], selection)
                        self.assertEqual(actual, expected)
                        skim = self.output / f"{system}_{definition}_{region}.root"
                        receipt = write_recoil_skim([retained], skim, selection)
                        self.assertEqual(receipt["selected_pairs"], expected["selected_pairs"])
                        with uproot.open(skim) as root:
                            rows = root["recoilPairs"].arrays(library="np")
                        self.assertTrue(np.all(rows["event_id_lo"] == 1002))
                        sumw, _ = np.histogram(rows["xjgamma"], bins=actual["axis"]["edges"], weights=rows["event_weight"])
                        sumw2, _ = np.histogram(rows["xjgamma"], bins=actual["axis"]["edges"], weights=rows["event_weight"] ** 2)
                        self.assertEqual(sumw.tolist(), expected["sumw"])
                        self.assertEqual(sumw2.tolist(), expected["sumw2"])
                for radius in (0.3, 0.4):
                    with self.subTest(system=system, definition=definition, radius=radius):
                        self.assertEqual(
                            purity_counts([retained], non_tight_definition=definition, isolation_radius=radius),
                            purity_counts([subset], non_tight_definition=definition, isolation_radius=radius),
                        )

            for dimension in ("1D", "2D"):
                with self.subTest(system=system, dimension=dimension):
                    config = ResponseBuildConfig(system=system, dimension=dimension)
                    actual = build_response([retained], config)
                    expected = build_response([subset], config)
                    for name, value in vars(actual).items():
                        if name == "provenance":
                            continue  # Different source inventories are intentional.
                        if isinstance(value, np.ndarray):
                            np.testing.assert_array_equal(value, getattr(expected, name))
                        elif isinstance(value, dict):
                            self.assertEqual(value.keys(), getattr(expected, name).keys())
                            for key, item in value.items():
                                np.testing.assert_array_equal(item, getattr(expected, name)[key])
                        else:
                            self.assertEqual(value, getattr(expected, name))

    def test_terminal_inventory_is_distinct_from_plotted_event_selection(self) -> None:
        specification, contract = self._specification("auau")
        for rows in specification["trees"].values():
            for row in rows:
                if row.get("event_id_lo") == 1001 and "terminal_status" in row:
                    row["terminal_status"] = 2
        source = self._write("mixed.root", "auau", specification, contract)
        selection = RecoilSelection(region="inclusive")
        histogram = self.output / "histogram.json"
        receipt_path = self.output / "histogram.receipt.json"
        dataset = ROOT / "config/datasets/auau_engineering_fixture.yaml"
        write_json(histogram, recoil_histogram([source], selection))
        receipt = write_histogram_receipt(
            input_paths=[source], histogram_path=histogram, selection=selection,
            dataset_manifest_path=dataset, receipt_path=receipt_path,
        )
        self.assertEqual(receipt["dataset_observation"]["event_count"], 2)
        self.assertEqual(receipt["dataset_observation"]["accepted_events"], 1)
        self.assertEqual(receipt["dataset_observation"]["retained_terminal_events_excluded"], 1)
        plot = compile_plot_contract(
            histogram_path=histogram, histogram_receipt_path=receipt_path,
            dataset_manifest_path=dataset,
        )
        self.assertIn("2 input events", plot["annotations"]["root_tlatex"]["dataset"])
        self.assertIn("Accepted producer events", plot["annotations"]["root_tlatex"]["cuts"])

    def test_all_terminal_input_has_empty_offline_products(self) -> None:
        for system in ("pp", "auau"):
            specification, contract = self._specification(system)
            for rows in specification["trees"].values():
                for row in rows:
                    if "terminal_status" in row:
                        row["terminal_status"] = 2
            source = self._write(f"{system}_all_terminal.root", system, specification, contract)
            for region in ("inclusive", "A"):
                selection = RecoilSelection(region=region)
                histogram = recoil_histogram([source], selection)
                self.assertEqual(histogram["selected_events"], 0)
                self.assertEqual(histogram["selected_pairs"], 0)
                self.assertEqual(sum(histogram["sumw"]), 0.0)
                skim = self.output / f"{system}_{region}_empty.root"
                write_recoil_skim([source], skim, selection)
                with uproot.open(skim) as root:
                    self.assertEqual(root["recoilPairs"].num_entries, 0)
            purity = purity_counts([source], non_tight_definition="bounded")
            self.assertEqual(sum(purity["weighted_counts"].values()), 0.0)
            for dimension in ("1D", "2D"):
                bundle = build_response([source], ResponseBuildConfig(system=system, dimension=dimension))
                self.assertEqual(float(bundle.reco.sum()), 0.0)
                self.assertEqual(float(bundle.truth.sum()), 0.0)
                self.assertEqual(bundle.provenance["observed"]["accepted_events"], 0)
                self.assertEqual(bundle.provenance["observed"]["retained_terminal_events_excluded"], 2)

    def test_truth_outside_classification_support_is_a_reco_boundary(self) -> None:
        specification, contract = self._specification("auau")
        specification["trees"]["truthPhotons"][0]["truth_photon_pt"] = 45.0
        source = self._write("boundary.root", "auau", specification, contract)
        bundle = build_response(
            [source], ResponseBuildConfig(system="auau", dimension="1D")
        )
        key = Category.PHOTON_HIGH_FEED_IN.value
        self.assertEqual(float(bundle.boundary_reco[key].sum()), 1.0)
        self.assertEqual(float(bundle.matrix.sum()), 1.0)
        self.assertEqual(float(bundle.fakes.sum()), 2.0)
        self.assertEqual(float(bundle.reco.sum()), 4.0)
        self.assertEqual(float(bundle.truth.sum()), 1.0)
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_missing_photon_truth_link_populates_unmatched_reco(self) -> None:
        specification, contract = self._specification("auau")
        photon = specification["trees"]["photons"][0]
        candidate_id = (photon["candidate_id_hi"], photon["candidate_id_lo"])
        photon["truth_matched"] = 0
        photon["truth_barcode"] = -1
        photon["truth_generator_occurrence_embedding_id"] = -1
        for pair in specification["trees"]["photonJets"]:
            if (pair["candidate_id_hi"], pair["candidate_id_lo"]) == candidate_id:
                pair["truth_matched"] = 0
                pair["truth_barcode"] = -1
                pair["truth_generator_occurrence_embedding_id"] = -1
        specification["trees"]["recoTruthLinks"] = [
            row
            for row in specification["trees"]["recoTruthLinks"]
            if not (
                int(row["reco_type"]) == 1
                and (row["reco_id_hi"], row["reco_id_lo"]) == candidate_id
            )
        ]
        source = self._write("unmatched_reco.root", "auau", specification, contract)
        bundle = build_response(
            [source],
            ResponseBuildConfig(system="auau", dimension="1D"),
        )
        self.assertEqual(
            float(bundle.fake_causes[Category.UNMATCHED_RECO.value].sum()),
            2.0,
        )
        self.assertEqual(float(bundle.fakes.sum()), 3.0)
        self.assertEqual(float(bundle.misses.sum()), 1.0)
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_explicit_reco_fake_and_truth_miss_links_preserve_partition(self) -> None:
        specification, contract = self._specification("auau")
        photon = specification["trees"]["photons"][0]
        candidate_id = (photon["candidate_id_hi"], photon["candidate_id_lo"])
        photon["truth_matched"] = 0
        photon["truth_barcode"] = -1
        photon["truth_generator_occurrence_embedding_id"] = -1
        for pair in specification["trees"]["photonJets"]:
            if (pair["candidate_id_hi"], pair["candidate_id_lo"]) == candidate_id:
                pair["truth_matched"] = 0
                pair["truth_barcode"] = -1
                pair["truth_generator_occurrence_embedding_id"] = -1

        links = specification["trees"]["recoTruthLinks"]
        link_index = next(
            index
            for index, row in enumerate(links)
            if int(row["reco_type"]) == 1
            and (row["reco_id_hi"], row["reco_id_lo"]) == candidate_id
        )
        original = deepcopy(links[link_index])
        reco_fake = deepcopy(original)
        reco_fake["truth_type"] = 0
        reco_fake["truth_id_hi"] = 0
        reco_fake["truth_id_lo"] = 0
        reco_fake["truth_index"] = -1
        reco_fake["link_class"] = 1
        truth_miss = deepcopy(original)
        truth_miss["link_id_lo"] = 987654320
        truth_miss["reco_type"] = 0
        truth_miss["reco_id_hi"] = 0
        truth_miss["reco_id_lo"] = 0
        truth_miss["reco_index"] = -1
        truth_miss["link_class"] = 2
        links[link_index] = reco_fake
        links.append(truth_miss)

        source = self._write("explicit_fake_miss.root", "auau", specification, contract)
        bundle = build_response(
            [source], ResponseBuildConfig(system="auau", dimension="1D")
        )
        self.assertEqual(
            float(bundle.fake_causes[Category.UNMATCHED_RECO.value].sum()),
            2.0,
        )
        self.assertEqual(float(bundle.misses.sum()), 1.0)
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_source_preserved_photon_many_to_one_and_jet_candidate_links(self) -> None:
        specification, contract = self._specification("auau")
        photons = specification["trees"]["photons"]
        first, second = photons[0], photons[1]
        second["truth_matched"] = 1
        second["truth_barcode"] = first["truth_barcode"]
        second["truth_generator_occurrence_embedding_id"] = first[
            "truth_generator_occurrence_embedding_id"
        ]
        second_id = (second["candidate_id_hi"], second["candidate_id_lo"])
        for pair in specification["trees"]["photonJets"]:
            if (pair["candidate_id_hi"], pair["candidate_id_lo"]) == second_id:
                pair["truth_matched"] = 1
                pair["truth_barcode"] = second["truth_barcode"]
                pair["truth_generator_occurrence_embedding_id"] = second[
                    "truth_generator_occurrence_embedding_id"
                ]

        links = specification["trees"]["recoTruthLinks"]
        photon_match = deepcopy(
            next(
                row
                for row in links
                if int(row["reco_type"]) == 1 and row["event_id_lo"] == 1001
            )
        )
        photon_match["link_id_lo"] = 987654318
        photon_match["reco_id_hi"] = second["candidate_id_hi"]
        photon_match["reco_id_lo"] = second["candidate_id_lo"]
        photon_match["reco_index"] = 1
        links.append(photon_match)

        jet_candidate = deepcopy(
            next(
                row
                for row in links
                if int(row["reco_type"]) == 2 and row["event_id_lo"] == 1001
            )
        )
        jet_candidate["link_id_lo"] = 987654319
        jet_candidate["link_class"] = 5
        links.append(jet_candidate)

        source = self._write("source_link_topologies.root", "auau", specification, contract)
        bundle = build_response(
            [source], ResponseBuildConfig(system="auau", dimension="1D")
        )
        self.assertEqual(float(bundle.matrix.sum()), 2.0)
        self.assertEqual(float(bundle.fakes.sum()), 2.0)
        self.assertTrue(bundle.conservation()["closes_rtol_1e-10_atol_1e-10"])

    def test_ambiguous_truth_signal_fails_closed(self) -> None:
        specification, contract = self._specification("auau")
        overlay = specification["trees"]["truthPhotons"][-1]
        overlay["source_role"] = 1
        overlay["truth_photon_pt"] = 20.0
        overlay["truth_isolation"] = 0.5
        source = self._write("ambiguous_truth.root", "auau", specification, contract)
        with self.assertRaisesRegex(ValueError, "ambiguous truth-signal photons"):
            build_response([source], ResponseBuildConfig(system="auau"))

    def test_duplicate_typed_link_fails_closed(self) -> None:
        specification, contract = self._specification("auau")
        duplicate = deepcopy(specification["trees"]["recoTruthLinks"][0])
        duplicate["link_id_lo"] = 987654321
        specification["trees"]["recoTruthLinks"].append(duplicate)
        source = self._write("duplicate_link.root", "auau", specification, contract)
        with self.assertRaisesRegex(ValueError, "ambiguous typed reco link"):
            build_response([source], ResponseBuildConfig(system="auau"))

    def test_pair_identity_cannot_disagree_with_local_object_indices(self) -> None:
        specification, contract = self._specification("auau")
        specification["trees"]["photonJets"][0]["jet_id_lo"] = 3103
        source = self._write("cross_event_pair_jet.root", "auau", specification, contract)
        with self.assertRaisesRegex(
            ValueError,
            "photonJets jet field disagrees with its local index|"
            "pair jet index and identity disagree",
        ):
            build_response([source], ResponseBuildConfig(system="auau"))

    def test_builder_rechecks_pair_witnesses_against_normalized_objects(self) -> None:
        specification, contract = self._specification("auau")
        specification["trees"]["photonJets"][0]["jet_eta"] = 99.0
        source = self._write("wrong_pair_copy.root", "auau", specification, contract)
        with mock.patch.object(response_builder_module, "validate_trees", return_value=[]):
            with self.assertRaisesRegex(
                ValueError,
                "photonJets jet_eta disagrees with normalized objects",
            ):
                build_response([source], ResponseBuildConfig(system="auau"))

    def test_builder_rejects_a_missing_photon_jet_combination(self) -> None:
        specification, contract = self._specification("auau")
        specification["trees"]["photonJets"].pop(0)
        source = self._write("missing_pair.root", "auau", specification, contract)
        with mock.patch.object(response_builder_module, "validate_trees", return_value=[]):
            with self.assertRaisesRegex(
                ValueError,
                "does not exactly cover the photon/jet Cartesian product",
            ):
                build_response([source], ResponseBuildConfig(system="auau"))

    def test_unproved_wrong_object_witness_fails_closed(self) -> None:
        specification, contract = self._specification("auau")
        specification["trees"]["photonJets"][0]["wrong_recoil_class"] = 1
        source = self._write("wrong_recoil.root", "auau", specification, contract)
        with self.assertRaisesRegex(ValueError, "golden response semantics"):
            build_response([source], ResponseBuildConfig(system="auau"))

    def test_wrong_object_witness_fails_even_when_no_region_leader_exists(self) -> None:
        specification, contract = self._specification("pp")
        specification["trees"]["photonJets"][0]["wrong_recoil_class"] = 1
        source = self._write("no_leader_wrong_recoil.root", "pp", specification, contract)
        with self.assertRaisesRegex(ValueError, "golden response semantics"):
            build_response([source], ResponseBuildConfig(system="pp"))

    def test_system_and_duplicate_input_ambiguity_fail_closed(self) -> None:
        pp = FIXTURES / "photonjet_trees_pp.root"
        with self.assertRaisesRegex(ValueError, "centrality"):
            build_response([pp], ResponseBuildConfig(system="auau"))
        with self.assertRaisesRegex(ValueError, "duplicate response input path"):
            build_response([pp, pp], ResponseBuildConfig(system="pp"))

    def test_event_source_identity_cannot_repeat_across_distinct_inputs(self) -> None:
        specification, contract = self._specification("auau")
        for rows in specification["trees"].values():
            for row in rows:
                if row.get("event_id_lo") == 1001 and "event_weight" in row:
                    row["event_weight"] = 1.5
        second = self._write("second_shard.root", "auau", specification, contract)
        with self.assertRaisesRegex(
            ValueError,
            "identity repeats across parts|event/source identity is repeated",
        ):
            build_response(
                [FIXTURES / "photonjet_trees_auau.root", second],
                ResponseBuildConfig(system="auau"),
            )

    def test_object_identity_cannot_repeat_across_distinct_input_parts(self) -> None:
        specification, contract = self._specification("auau")
        for rows in specification["trees"].values():
            for row in rows:
                if "source_file_index" in row:
                    row["source_file_index"] = int(row["source_file_index"]) + 1
                if "event_id_lo" in row:
                    row["event_id_lo"] = int(row["event_id_lo"]) + 10000
        second = self._write("repeated_objects.root", "auau", specification, contract)
        with self.assertRaisesRegex(ValueError, "identity repeats across parts"):
            build_response(
                [FIXTURES / "photonjet_trees_auau.root", second],
                ResponseBuildConfig(system="auau"),
            )

    def test_artifact_replay_is_byte_exact_and_path_portable(self) -> None:
        source = FIXTURES / "photonjet_trees_auau.root"
        stem = self.output / "response"
        config = ResponseBuildConfig(system="auau", dimension="2D")
        first = write_response_artifacts([source], stem, config)
        first_bytes = {
            suffix: stem.with_suffix(suffix).read_bytes()
            for suffix in (".npz", ".json", ".receipt.json")
        }
        second = write_response_artifacts([source], stem, config)
        second_bytes = {
            suffix: stem.with_suffix(suffix).read_bytes()
            for suffix in (".npz", ".json", ".receipt.json")
        }
        self.assertEqual(first, second)
        self.assertEqual(first_bytes, second_bytes)
        self.assertTrue(first["diagnostics"]["conservation"]["closes_rtol_1e-10_atol_1e-10"])
        serialized = json.dumps(first, sort_keys=True)
        self.assertNotIn(self.temporary.name, serialized)
        self.assertEqual(first["unresolved"]["golden_mini_dst_differential"], "NOT_RUN")

    def test_custom_receipt_directory_is_supported(self) -> None:
        source = FIXTURES / "photonjet_trees_auau.root"
        stem = self.output / "payloads" / "response"
        receipt = self.output / "receipts" / "accepted.json"
        payload = write_response_artifacts(
            [source],
            stem,
            ResponseBuildConfig(system="auau", dimension="2D"),
            receipt_path=receipt,
        )
        self.assertTrue(receipt.is_file())
        self.assertEqual(json.loads(receipt.read_text(encoding="utf-8")), payload)
        self.assertEqual(
            [row["name"] for row in payload["outputs"]],
            ["response.npz", "response.json"],
        )

    def test_interrupted_republication_cannot_leave_a_stale_receipt(self) -> None:
        source = FIXTURES / "photonjet_trees_auau.root"
        stem = self.output / "response"
        published_stem = stem.resolve()
        config = ResponseBuildConfig(system="auau", dimension="2D")
        write_response_artifacts([source], stem, config)
        receipt = published_stem.with_suffix(".receipt.json")
        self.assertTrue(receipt.is_file())

        real_replace = __import__("os").replace
        payload_replacements = 0

        def fail_during_payload_publication(source_path: Path, destination_path: Path) -> None:
            nonlocal payload_replacements
            destination = Path(destination_path)
            if destination in {
                published_stem.with_suffix(".npz"),
                published_stem.with_suffix(".json"),
            }:
                payload_replacements += 1
                if payload_replacements == 2:
                    raise OSError("injected publication interruption")
            real_replace(source_path, destination_path)

        with mock.patch(
            "photonjet.analysis.response_builder.os.replace",
            side_effect=fail_during_payload_publication,
        ):
            with self.assertRaisesRegex(OSError, "injected publication interruption"):
                write_response_artifacts([source], stem, config)
        self.assertFalse(receipt.exists())

    def test_input_change_during_read_fails_before_receipted_result(self) -> None:
        specification, contract = self._specification("auau")
        source = self._write("mutable.root", "auau", specification, contract)
        for row in specification["trees"]["events"]:
            row["event_weight"] = float(row["event_weight"]) + 1.0
        for tree_name in ("eventTree", "photons", "jets", "photonJets"):
            for row in specification["trees"][tree_name]:
                row["event_weight"] = float(row["event_weight"]) + 1.0
        replacement = self._write("replacement.root", "auau", specification, contract)
        real_validate = response_builder_module.validate_trees

        def validate_then_replace(*args, **kwargs):
            result = real_validate(*args, **kwargs)
            shutil.copyfile(replacement, source)
            return result

        with mock.patch.object(
            response_builder_module,
            "validate_trees",
            side_effect=validate_then_replace,
        ):
            with self.assertRaisesRegex(
                ValueError,
                "response input changed while it was being validated or read",
            ):
                build_response([source], ResponseBuildConfig(system="auau"))

    def test_concurrent_writer_for_the_same_output_fails_closed(self) -> None:
        source = FIXTURES / "photonjet_trees_auau.root"
        stem = (self.output / "response").resolve()
        lock_path = stem.with_suffix(".npz").with_name(".response.npz.lock")
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        with lock_path.open("a+b") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(
                RuntimeError,
                "another response writer owns output path",
            ):
                write_response_artifacts(
                    [source],
                    stem,
                    ResponseBuildConfig(system="auau"),
                )
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def test_cli_response_build_uses_the_same_stage(self) -> None:
        stem = self.output / "cli_response"
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            status = cli_main(
                [
                    "response",
                    "build",
                    "--input",
                    str(FIXTURES / "photonjet_trees_auau.root"),
                    "--system",
                    "auau",
                    "--dimension",
                    "2D",
                    "--region",
                    "A",
                    "--output-stem",
                    str(stem),
                ]
            )
        self.assertEqual(status, 0)
        summary = json.loads(stdout.getvalue())
        self.assertEqual(
            summary["status"],
            "TREE_CONSTRUCTION_COMPLETE__GOLDEN_EQUIVALENCE_NOT_RUN",
        )
        self.assertEqual(summary["outputs"], ["cli_response.npz", "cli_response.json"])
        self.assertTrue(stem.with_suffix(".receipt.json").is_file())


if __name__ == "__main__":
    unittest.main()
