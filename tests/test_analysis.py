from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection, recoil_histogram, write_recoil_skim


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


class AnalysisTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = FIXTURES / "photonjet_trees_pp.root"

    def _tampered_photons(
        self,
        branch: str,
        value: int | float,
        *,
        row: int = 0,
        source: Path | None = None,
        dtype: str | None = None,
        additional_updates: dict[str, int | float] | None = None,
    ) -> Path:
        input_path = self.output if source is None else source
        with uproot.open(input_path) as root:
            arrays = {
                name: np.array(values, copy=True)
                for name, values in root["photons"].arrays(library="np").items()
            }
        if dtype is not None:
            arrays[branch] = arrays[branch].astype(dtype)
        arrays[branch][row] = value
        for name, updated in (additional_updates or {}).items():
            arrays[name][row] = updated
        output = Path(self.temporary.name) / f"tampered_{branch}.root"
        with uproot.recreate(output) as root:
            root["photons"] = arrays
        return output

    def _reordered_photons(self, order: list[int]) -> Path:
        with uproot.open(self.output) as root:
            arrays = {
                name: np.asarray(values)[order]
                for name, values in root["photons"].arrays(library="np").items()
            }
        output = Path(self.temporary.name) / "reordered_photons.root"
        with uproot.recreate(output) as root:
            root["photons"] = arrays
        return output

    def _photon_subset(self, rows: list[int], name: str) -> Path:
        with uproot.open(self.output) as root:
            arrays = {
                branch: np.asarray(values)[rows]
                for branch, values in root["photons"].arrays(library="np").items()
            }
        output = Path(self.temporary.name) / name
        with uproot.recreate(output) as root:
            root["photons"] = arrays
        return output

    def test_histogram_is_deterministic_and_complete(self) -> None:
        selection = RecoilSelection(region="inclusive")
        first = recoil_histogram([self.output], selection)
        second = recoil_histogram([self.output], selection)
        self.assertEqual(first, second)
        self.assertEqual(len(first["sumw"]), len(first["axis"]["edges"]) - 1)
        self.assertEqual(first["selected_pairs"], int(sum(first["sumw"])))

    def test_event_leading_skim_preserves_selection(self) -> None:
        selection = RecoilSelection(region="A")
        skim = Path(self.temporary.name) / "skim.root"
        receipt = write_recoil_skim([self.output], skim, selection)
        with uproot.open(skim) as root:
            self.assertEqual(root["recoilPairs"].num_entries, receipt["selected_pairs"])
            values = root["recoilPairs"].arrays(["photon_et", "jet_pt", "delta_phi"], library="np")
        self.assertTrue(np.all(values["photon_et"] >= selection.photon_et_min))
        self.assertTrue(np.all(values["photon_et"] < selection.photon_et_max))
        self.assertTrue(np.all(values["jet_pt"] > selection.jet_pt_min))
        self.assertTrue(np.all(values["delta_phi"] > selection.delta_phi_min))
        self.assertEqual(receipt["output"]["name"], "skim.root")
        self.assertEqual(receipt["output"]["role"], "recoil_skim")
        serialized = json.dumps(receipt, sort_keys=True)
        self.assertNotIn(self.temporary.name, serialized)
        root_bytes = skim.read_bytes()
        for marker in (
            b"/" + b"Users/",
            b"/" + b"home/",
            b"/private" + b"/" + b"tmp/",
            b"/" + b"tmp/",
            b"/" + b"sphenix/",
            b"/" + b"gpfs",
        ):
            self.assertNotIn(marker, root_bytes)

    def test_purity_boundaries_are_explicit(self) -> None:
        bounded = purity_counts(
            [self.output],
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        complement = purity_counts(
            [self.output],
            non_tight_definition="complement",
            isolation_radius=0.4,
        )
        self.assertEqual(set(bounded["weighted_counts"]), set("ABCD"))
        self.assertGreaterEqual(
            complement["weighted_counts"]["C"],
            bounded["weighted_counts"]["C"],
        )

    def test_purity_rejects_mixed_terminal_state_within_one_event(self) -> None:
        tampered = self._tampered_photons("terminal_status", 2)
        with self.assertRaisesRegex(ValueError, "different terminal status"):
            purity_counts([tampered], non_tight_definition="complement")

    def test_purity_rejects_fractional_terminal_status(self) -> None:
        tampered = self._tampered_photons("terminal_status", 0.2, dtype="float64")
        with self.assertRaisesRegex(ValueError, "non-negative int32"):
            purity_counts([tampered], non_tight_definition="complement")

    def test_purity_fixture_payloads_are_sealed(self) -> None:
        expected = {
            ("pp", "bounded"): [0, 0, 0, 0],
            ("pp", "complement"): [0, 0, 2, 2],
            ("auau", "bounded"): [2, 2, 0, 0],
            ("auau", "complement"): [2, 2, 0, 0],
        }
        for (system, definition), values in expected.items():
            for radius in (0.3, 0.4):
                with self.subTest(
                    system=system,
                    definition=definition,
                    radius=radius,
                ):
                    payload = purity_counts(
                        [FIXTURES / f"photonjet_trees_{system}.root"],
                        non_tight_definition=definition,
                        isolation_radius=radius,
                    )
                    self.assertEqual(
                        list(payload["event_counts"].values()),
                        values,
                    )
                    self.assertEqual(
                        list(payload["weighted_counts"].values()),
                        [float(value) for value in values],
                    )

    def test_purity_rejects_nonbinary_classification_witnesses(self) -> None:
        tampered = self._tampered_photons("bdt_is_not_tight", 2)
        with self.assertRaisesRegex(ValueError, "complement BDT witness is not binary"):
            purity_counts(
                [tampered],
                non_tight_definition="complement",
                isolation_radius=0.4,
            )

    def test_purity_rejects_nonbinary_bounded_witnesses(self) -> None:
        tampered = self._tampered_photons("bdt_is_nontight", 2)
        with self.assertRaisesRegex(ValueError, "bounded BDT witness is not binary"):
            purity_counts(
                [tampered],
                non_tight_definition="bounded",
                isolation_radius=0.4,
            )

    def test_purity_ignores_inactive_non_tight_witnesses(self) -> None:
        bounded_nominal = purity_counts(
            [self.output],
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        bounded_tampered = purity_counts(
            [self._tampered_photons("bdt_is_not_tight", 2)],
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        self.assertEqual(bounded_tampered, bounded_nominal)

        complement_nominal = purity_counts(
            [self.output],
            non_tight_definition="complement",
            isolation_radius=0.4,
        )
        complement_tampered = purity_counts(
            [self._tampered_photons("bdt_is_nontight", 2)],
            non_tight_definition="complement",
            isolation_radius=0.4,
        )
        self.assertEqual(complement_tampered, complement_nominal)

    def test_purity_rejects_nonbinary_isolation_witnesses(self) -> None:
        tampered = self._tampered_photons("iso_r03_pass", -7)
        with self.assertRaisesRegex(ValueError, "isolation r03 witness is not binary"):
            purity_counts(
                [tampered],
                non_tight_definition="bounded",
                isolation_radius=0.3,
            )

    def test_purity_rejects_binary_witnesses_that_disagree_with_values(self) -> None:
        cases = (
            (
                "bdt_is_tight",
                1,
                "bounded",
                "tight BDT witness differs",
            ),
            (
                "bdt_is_not_tight",
                0,
                "complement",
                "complement BDT witness differs",
            ),
            (
                "iso_r04_pass",
                0,
                "bounded",
                "isolation witness differs",
            ),
        )
        for branch, value, definition, message in cases:
            with self.subTest(branch=branch):
                tampered = self._tampered_photons(branch, value)
                with self.assertRaisesRegex(ValueError, message):
                    purity_counts(
                        [tampered],
                        non_tight_definition=definition,
                        isolation_radius=0.4,
                    )

    def test_purity_rejects_duplicate_event_ordinals(self) -> None:
        tampered = self._tampered_photons(
            "photon_encounter_ordinal",
            0,
            row=1,
        )
        with self.assertRaisesRegex(ValueError, "encounter ordinal is duplicated"):
            purity_counts(
                [tampered],
                non_tight_definition="complement",
                isolation_radius=0.4,
            )

    def test_purity_rejects_fractional_event_ordinals(self) -> None:
        tampered = self._tampered_photons(
            "photon_encounter_ordinal",
            0.5,
            dtype="float64",
        )
        with self.assertRaisesRegex(ValueError, "encounter ordinal is not an integer"):
            purity_counts(
                [tampered],
                non_tight_definition="complement",
                isolation_radius=0.4,
            )

    def test_purity_rejects_an_event_split_across_input_files(self) -> None:
        first = self._photon_subset([0], "split_first.root")
        second = self._photon_subset([1], "split_second.root")
        with self.assertRaisesRegex(ValueError, "event are not contiguous"):
            purity_counts(
                [first, second],
                non_tight_definition="complement",
                isolation_radius=0.4,
            )

    def test_purity_rejects_noncontiguous_event_rows(self) -> None:
        tampered = self._reordered_photons([0, 2, 1, 3])
        with self.assertRaisesRegex(ValueError, "event are not contiguous"):
            purity_counts(
                [tampered],
                non_tight_definition="complement",
                isolation_radius=0.4,
            )

    def test_purity_leader_eligibility_rejects_out_of_acceptance_photons(self) -> None:
        source = FIXTURES / "photonjet_trees_auau.root"
        for branch, value in (("photon_et", 36.0), ("photon_eta", 0.8)):
            with self.subTest(branch=branch):
                tampered = self._tampered_photons(
                    branch,
                    value,
                    row=0,
                    source=source,
                )
                payload = purity_counts(
                    [tampered],
                    non_tight_definition="bounded",
                    isolation_radius=0.4,
                )
                self.assertEqual(payload["event_counts"]["A"], 1)
                self.assertEqual(payload["event_counts"]["B"], 2)

    def test_purity_ignores_classification_witnesses_outside_acceptance(self) -> None:
        source = FIXTURES / "photonjet_trees_auau.root"
        tampered = self._tampered_photons(
            "photon_et",
            36.0,
            row=0,
            source=source,
            additional_updates={"bdt_is_tight": 2, "iso_r04_pass": -7},
        )
        payload = purity_counts(
            [tampered],
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        self.assertEqual(payload["event_counts"]["A"], 1)


if __name__ == "__main__":
    unittest.main()
