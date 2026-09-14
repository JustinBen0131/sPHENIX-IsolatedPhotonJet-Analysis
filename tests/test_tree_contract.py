from __future__ import annotations

from copy import deepcopy
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

from photonjet.io.tree_validation import validate


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
FIXTURE_GENERATOR_SPEC = importlib.util.spec_from_file_location(
    "photonjet_engineering_fixture_generator",
    ROOT / "tools/generate_engineering_fixtures.py",
)
if FIXTURE_GENERATOR_SPEC is None or FIXTURE_GENERATOR_SPEC.loader is None:
    raise RuntimeError("engineering-fixture generator cannot be loaded")
FIXTURE_GENERATOR = importlib.util.module_from_spec(FIXTURE_GENERATOR_SPEC)
FIXTURE_GENERATOR_SPEC.loader.exec_module(FIXTURE_GENERATOR)


class TreeContractTest(unittest.TestCase):
    def test_tree_validator_import_is_order_independent(self) -> None:
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(ROOT / "FinalAnalysis")
        completed = subprocess.run(
            [
                sys.executable,
                "-c",
                "from photonjet.io.tree_validation import validate; "
                "from photonjet.analysis.response_builder import build_response",
            ],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_pp_tree_contract(self) -> None:
        with uproot.open(FIXTURES / "photonjet_trees_pp.root") as root:
            centrality = root["events"]["centrality"].array(library="np")
        self.assertTrue(np.all(centrality == -1.0))
        lines = validate(
            [FIXTURES / "photonjet_trees_pp.root"],
            model_input_count=11,
        )
        self.assertIn("tree_schema_exact=PASS", lines)
        self.assertIn("eventTree_flat_equivalence=PASS", lines)
        self.assertIn("truth_link_type_event_index=PASS", lines)
        self.assertIn("truth_occurrence_identity=PASS", lines)
        self.assertIn(
            "counts=events:2,eventTree:2,photons:4,jets:4,photonJets:8,"
            "truthPhotons:3,truthJets:2,recoTruthLinks:4",
            lines,
        )

    def test_fixture_exercises_occurrence_safe_duplicate_barcodes(self) -> None:
        with uproot.open(FIXTURES / "photonjet_trees_pp.root") as root:
            rows = root["truthPhotons"].arrays(
                [
                    "event_id_lo",
                    "generator_barcode",
                    "generator_occurrence_embedding_id",
                ],
                library="np",
            )
        event_mask = rows["event_id_lo"] == 1001
        self.assertEqual(rows["generator_barcode"][event_mask].tolist(), [11, 11])
        self.assertEqual(
            rows["generator_occurrence_embedding_id"][event_mask].tolist(),
            [2, 1],
        )

    def test_duplicate_occurrence_and_barcode_in_one_event_fails_closed(self) -> None:
        records = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(encoding="utf-8")
        )
        contract = FIXTURE_GENERATOR.load_contract()
        specification = deepcopy(
            FIXTURE_GENERATOR.materialize_records(records, contract)["systems"]["pp"]
        )
        specification["trees"]["truthPhotons"][-1][
            "generator_occurrence_embedding_id"
        ] = 2
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "duplicate.root"
            FIXTURE_GENERATOR.write_fixture(
                output=output,
                system="pp",
                specification=specification,
                contract=contract,
            )
            with self.assertRaisesRegex(
                ValueError,
                "duplicate generator occurrence/barcode identity",
            ):
                validate([output], model_input_count=11)

    def test_reco_occurrence_must_resolve_to_the_typed_truth_link(self) -> None:
        records = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(encoding="utf-8")
        )
        contract = FIXTURE_GENERATOR.load_contract()
        specification = deepcopy(
            FIXTURE_GENERATOR.materialize_records(records, contract)["systems"]["pp"]
        )
        specification["trees"]["photons"][0][
            "truth_generator_occurrence_embedding_id"
        ] = 1
        candidate_id = (
            specification["trees"]["photons"][0]["candidate_id_hi"],
            specification["trees"]["photons"][0]["candidate_id_lo"],
        )
        for pair in specification["trees"]["photonJets"]:
            if (pair["candidate_id_hi"], pair["candidate_id_lo"]) == candidate_id:
                pair["truth_generator_occurrence_embedding_id"] = 1
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "wrong_occurrence.root"
            FIXTURE_GENERATOR.write_fixture(
                output=output,
                system="pp",
                specification=specification,
                contract=contract,
            )
            with self.assertRaisesRegex(
                ValueError,
                "occurrence/barcode disagrees with its typed truth link",
            ):
                validate([output], model_input_count=11)

    def test_pair_object_identities_must_match_their_local_indices(self) -> None:
        records = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(encoding="utf-8")
        )
        contract = FIXTURE_GENERATOR.load_contract()
        for side, field, replacement in (
            ("photon", "candidate_id_lo", 2103),
            ("photon", "photon_et", 19.5),
            ("jet", "jet_id_lo", 3103),
            ("jet", "jet_eta", 99.0),
        ):
            with self.subTest(side=side), tempfile.TemporaryDirectory() as temporary:
                specification = deepcopy(
                    FIXTURE_GENERATOR.materialize_records(records, contract)["systems"][
                        "auau"
                    ]
                )
                specification["trees"]["photonJets"][0][field] = replacement
                output = Path(temporary) / f"wrong_{side}_identity.root"
                FIXTURE_GENERATOR.write_fixture(
                    output=output,
                    system="auau",
                    specification=specification,
                    contract=contract,
                )
                with self.assertRaisesRegex(
                    ValueError,
                    rf"photonJets {side} field disagrees with its local index.*{field}",
                ):
                    validate([output], model_input_count=14)

    def test_flat_event_fields_must_match_their_parent_event(self) -> None:
        records = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(encoding="utf-8")
        )
        contract = FIXTURE_GENERATOR.load_contract()
        specification = deepcopy(
            FIXTURE_GENERATOR.materialize_records(records, contract)["systems"]["auau"]
        )
        specification["trees"]["photons"][0]["event_weight"] = 99.0
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "wrong_parent_event.root"
            FIXTURE_GENERATOR.write_fixture(
                output=output,
                system="auau",
                specification=specification,
                contract=contract,
            )
            with self.assertRaisesRegex(
                ValueError,
                "photons event field disagrees with its parent: event_weight",
            ):
                validate([output], model_input_count=14)

    def test_reco_truth_link_class_and_null_topology_fail_closed(self) -> None:
        records = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(encoding="utf-8")
        )
        contract = FIXTURE_GENERATOR.load_contract()
        for defect, mutate, message in (
            (
                "nonnul_absent_truth",
                lambda link: link.update(
                    truth_type=0,
                    truth_index=-1,
                    link_class=1,
                ),
                "absent truth target is not the canonical null",
            ),
            (
                "fake_with_truth",
                lambda link: link.update(link_class=1),
                "link_class topology differs",
            ),
        ):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as temporary:
                specification = deepcopy(
                    FIXTURE_GENERATOR.materialize_records(records, contract)["systems"][
                        "auau"
                    ]
                )
                mutate(specification["trees"]["recoTruthLinks"][0])
                output = Path(temporary) / f"{defect}.root"
                FIXTURE_GENERATOR.write_fixture(
                    output=output,
                    system="auau",
                    specification=specification,
                    contract=contract,
                )
                with self.assertRaisesRegex(ValueError, message):
                    validate([output], model_input_count=14)

    def test_pair_rows_must_cover_every_local_photon_jet_combination(self) -> None:
        records = json.loads(
            (FIXTURES / "engineering_fixture_seed_v1.json").read_text(encoding="utf-8")
        )
        contract = FIXTURE_GENERATOR.load_contract()
        specification = deepcopy(
            FIXTURE_GENERATOR.materialize_records(records, contract)["systems"]["auau"]
        )
        specification["trees"]["photonJets"].pop(0)
        event_tree = specification["trees"]["eventTree"][0]
        event_tree["npairs"] -= 1
        for field in (
            "pair_photon_index",
            "pair_jet_index",
            "pair_delta_phi",
            "pair_xjgamma",
            "pair_recoil_state",
        ):
            event_tree[field].pop(0)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "missing_pair.root"
            FIXTURE_GENERATOR.write_fixture(
                output=output,
                system="auau",
                specification=specification,
                contract=contract,
            )
            with self.assertRaisesRegex(
                ValueError,
                "photonJets does not cover the local photon/jet Cartesian product",
            ):
                validate([output], model_input_count=14)

    def test_auau_tree_contract(self) -> None:
        lines = validate(
            [FIXTURES / "photonjet_trees_auau.root"],
            model_input_count=14,
        )
        self.assertIn("eventTree_executable_leaders=PASS", lines)
        self.assertIn("numeric_contract=PASS", lines)


if __name__ == "__main__":
    unittest.main()
