from __future__ import annotations

import unittest

from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.selection import (
    Predicate,
    SelectionProgram,
    compile_abcd_selection,
    compile_recoil_selection,
)


class SelectionContractTest(unittest.TestCase):
    def test_boundaries_are_executed_by_the_serialized_predicates(self) -> None:
        program = compile_recoil_selection(RecoilSelection(region="inclusive"))
        self.assertTrue(program.accepts("photon", {"photon_et": 15.0, "photon_eta": 0.0}))
        self.assertFalse(program.accepts("photon", {"photon_et": 35.0, "photon_eta": 0.0}))
        self.assertFalse(program.accepts("photon", {"photon_et": 20.0, "photon_eta": 0.7}))
        self.assertFalse(
            program.accepts(
                "recoil",
                {"jet_pt": 5.0, "jet_eta": 0.0, "jet_radius": 0.4, "delta_phi": 3.0},
            )
        )
        self.assertTrue(
            program.accepts(
                "recoil",
                {"jet_pt": 5.1, "jet_eta": 0.1, "jet_radius": 0.4, "delta_phi": 3.0},
            )
        )

    def test_irrelevant_non_tight_definition_is_pruned(self) -> None:
        region_a = compile_recoil_selection(
            RecoilSelection(region="A", non_tight_definition="bounded")
        )
        region_c = compile_recoil_selection(
            RecoilSelection(region="C", non_tight_definition="complement")
        )
        a_keys = {fact.key for fact in region_a.facts}
        c_values = {fact.key: fact.value for fact in region_c.facts}
        self.assertNotIn("photon.non_tight_definition", a_keys)
        self.assertEqual(c_values["photon.non_tight_definition"], "complement")
        self.assertEqual(region_c.leader_branch, "leader_C_r04_complement_index")

    def test_recoil_and_purity_compile_the_same_abcd_program(self) -> None:
        for region in "ABCD":
            for definition in ("bounded", "complement"):
                with self.subTest(region=region, definition=definition):
                    standalone = compile_abcd_selection(
                        region=region,
                        non_tight_definition=definition,
                        isolation_radius=0.4,
                    )
                    recoil = compile_recoil_selection(
                        RecoilSelection(
                            region=region,
                            non_tight_definition=definition,
                        )
                    )
                    recoil_photon_program = tuple(
                        predicate
                        for predicate in recoil.predicates
                        if predicate.stage in {"event", "photon", "photon_class"}
                    )
                    self.assertEqual(recoil_photon_program, standalone.predicates)
                    self.assertEqual(recoil.facts, standalone.facts)
                    self.assertEqual(recoil.leader_branch, standalone.leader_branch)

    def test_accepted_event_predicate_cannot_be_removed_or_redefined(self) -> None:
        for region in ("inclusive", "A", "B", "C", "D"):
            program = compile_recoil_selection(RecoilSelection(region=region))
            self.assertTrue(program.accepts_event(0))
            for status in (1, 2, 3, 4, 5):
                self.assertFalse(program.accepts_event(status))
            removed = program.to_dict()
            removed["predicates"] = [
                row for row in removed["predicates"] if row["stage"] != "event"
            ]
            with self.assertRaisesRegex(ValueError, "accepted-event predicate"):
                SelectionProgram.from_dict(removed)
            changed = program.to_dict()
            changed["predicates"][0]["value"] = 1.0
            with self.assertRaisesRegex(ValueError, "accepted-event predicate"):
                SelectionProgram.from_dict(changed)

    def test_terminal_status_is_not_coerced_to_success(self) -> None:
        program = compile_recoil_selection(RecoilSelection())
        for status in (False, True, 0.0, 0.2, "0", None, -1, 2**31, float("nan")):
            with self.subTest(status=status):
                with self.assertRaisesRegex(ValueError, "non-negative int32"):
                    program.accepts_event(status)

    def test_abcd_program_identity_tracks_every_relevant_choice(self) -> None:
        nominal = compile_abcd_selection(
            region="C",
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        variants = (
            compile_abcd_selection(
                region="D",
                non_tight_definition="bounded",
                isolation_radius=0.4,
            ),
            compile_abcd_selection(
                region="C",
                non_tight_definition="complement",
                isolation_radius=0.4,
            ),
            compile_abcd_selection(
                region="C",
                non_tight_definition="bounded",
                isolation_radius=0.3,
            ),
        )
        self.assertEqual(len({nominal.sha256, *(item.sha256 for item in variants)}), 4)

        tight_bounded = compile_abcd_selection(
            region="A",
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        tight_complement = compile_abcd_selection(
            region="A",
            non_tight_definition="complement",
            isolation_radius=0.4,
        )
        self.assertEqual(tight_bounded.sha256, tight_complement.sha256)

    def test_abcd_boundaries_are_strict_and_exhaustive(self) -> None:
        base = {
            "photon_bdt_score": 0.8,
            "photon_bdt_tight_threshold": 0.8,
            "photon_bdt_nontight_low_threshold": 0.5,
            "photon_bdt_nontight_high_threshold": 0.7,
            "photon_iso_r04": 4.0,
            "photon_iso_r04_threshold": 4.0,
            "photon_iso_r04_nonisolated_threshold": 7.0,
        }
        programs = {
            region: compile_abcd_selection(
                region=region,
                non_tight_definition="bounded",
                isolation_radius=0.4,
            )
            for region in "ABCD"
        }
        self.assertFalse(programs["A"].accepts("photon_class", base))

        tight_isolated = {**base, "photon_bdt_score": 0.800001, "photon_iso_r04": 3.999999}
        self.assertTrue(programs["A"].accepts("photon_class", tight_isolated))

        tight_nonisolated_boundary = {
            **base,
            "photon_bdt_score": 0.800001,
            "photon_iso_r04": 7.0,
        }
        self.assertFalse(
            programs["B"].accepts("photon_class", tight_nonisolated_boundary)
        )
        self.assertTrue(
            programs["B"].accepts(
                "photon_class",
                {**tight_nonisolated_boundary, "photon_iso_r04": 7.000001},
            )
        )

        for boundary in (0.5, 0.7):
            self.assertFalse(
                programs["C"].accepts(
                    "photon_class",
                    {**base, "photon_bdt_score": boundary, "photon_iso_r04": 1.0},
                )
            )
        self.assertTrue(
            programs["C"].accepts(
                "photon_class",
                {**base, "photon_bdt_score": 0.6, "photon_iso_r04": 1.0},
            )
        )

        complement = compile_abcd_selection(
            region="C",
            non_tight_definition="complement",
            isolation_radius=0.4,
        )
        self.assertTrue(
            complement.accepts(
                "photon_class",
                {**base, "photon_bdt_score": 0.8, "photon_iso_r04": 1.0},
            )
        )
        self.assertFalse(
            complement.accepts(
                "photon_class",
                {**base, "photon_bdt_score": 0.800001, "photon_iso_r04": 1.0},
            )
        )

    def test_r03_program_executes_leader_without_claiming_a_recorded_branch(self) -> None:
        program = compile_abcd_selection(
            region="A",
            non_tight_definition="bounded",
            isolation_radius=0.3,
        )
        records = [
            {
                "photon_et": 19.0,
                "photon_eta": 0.1,
                "photon_encounter_ordinal": 0,
                "photon_bdt_score": 0.90,
                "photon_bdt_tight_threshold": 0.80,
                "photon_iso_r03": 1.0,
                "photon_iso_r03_threshold": 4.0,
            },
            {
                "photon_et": 24.0,
                "photon_eta": 0.2,
                "photon_encounter_ordinal": 1,
                "photon_bdt_score": 0.91,
                "photon_bdt_tight_threshold": 0.80,
                "photon_iso_r03": 2.0,
                "photon_iso_r03_threshold": 4.0,
            },
        ]
        self.assertIsNone(program.leader_branch)
        self.assertEqual(program.choose_leader(records), 1)

    def test_serialized_program_rejects_duplicate_semantic_keys(self) -> None:
        value = compile_abcd_selection(
            region="A",
            non_tight_definition="bounded",
            isolation_radius=0.4,
        ).to_dict()
        value["facts"].append(dict(value["facts"][0]))
        with self.assertRaisesRegex(ValueError, "duplicate selection fact key"):
            SelectionProgram.from_dict(value)

    def test_serialized_program_rejects_fact_predicate_drift(self) -> None:
        value = compile_abcd_selection(
            region="A",
            non_tight_definition="bounded",
            isolation_radius=0.4,
        ).to_dict()
        for fact in value["facts"]:
            if fact["key"] == "photon.abcd_region":
                fact["value"] = "D"
            elif fact["key"] == "photon.isolation_state":
                fact["value"] = "nonisolated"
        with self.assertRaisesRegex(
            ValueError,
            "facts, predicates, and leader branch disagree",
        ):
            SelectionProgram.from_dict(value)

    def test_serialized_program_rejects_fictional_leader_rule(self) -> None:
        value = compile_abcd_selection(
            region="A",
            non_tight_definition="bounded",
            isolation_radius=0.4,
        ).to_dict()
        for fact in value["facts"]:
            if fact["key"] == "photon.leader_rule":
                fact["value"] = "highest_ordinal_then_lowest_et"
        with self.assertRaisesRegex(
            ValueError,
            "facts, predicates, and leader branch disagree",
        ):
            SelectionProgram.from_dict(value)

    def test_predicate_rejects_parameters_that_execution_ignores(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-executable parameters"):
            Predicate.from_dict(
                {
                    "key": "photon.bdt_class",
                    "stage": "photon_class",
                    "field": "photon_bdt_score",
                    "operator": "gt_field",
                    "reference": "photon_bdt_tight_threshold",
                    "value": 999.0,
                }
            )

    def test_every_cut_change_changes_the_program_identity(self) -> None:
        nominal = compile_recoil_selection(RecoilSelection())
        changed = compile_recoil_selection(RecoilSelection(jet_pt_min=6.0))
        self.assertNotEqual(nominal.sha256, changed.sha256)
        self.assertNotEqual(nominal.to_dict(), changed.to_dict())

    def test_abcd_predicates_execute_the_same_leader_state_that_is_serialized(self) -> None:
        program = compile_recoil_selection(RecoilSelection(region="A"))
        records = [
            {
                "photon_et": 19.0,
                "photon_eta": 0.1,
                "photon_encounter_ordinal": 0,
                "photon_bdt_score": 0.90,
                "photon_bdt_tight_threshold": 0.80,
                "photon_bdt_nontight_low_threshold": 0.50,
                "photon_bdt_nontight_high_threshold": 0.70,
                "photon_iso_r04": 1.0,
                "photon_iso_r04_threshold": 4.0,
                "photon_iso_r04_nonisolated_threshold": 7.0,
            },
            {
                "photon_et": 24.0,
                "photon_eta": 0.2,
                "photon_encounter_ordinal": 1,
                "photon_bdt_score": 0.91,
                "photon_bdt_tight_threshold": 0.80,
                "photon_bdt_nontight_low_threshold": 0.50,
                "photon_bdt_nontight_high_threshold": 0.70,
                "photon_iso_r04": 2.0,
                "photon_iso_r04_threshold": 4.0,
                "photon_iso_r04_nonisolated_threshold": 7.0,
            },
        ]
        self.assertEqual(program.choose_leader(records), 1)
        records[1]["photon_iso_r04"] = 8.0
        self.assertEqual(program.choose_leader(records), 0)

    def test_recoil_leader_eligibility_includes_photon_kinematics(self) -> None:
        program = compile_recoil_selection(RecoilSelection(region="A"))
        records = [
            {
                "photon_et": 36.0,
                "photon_eta": 0.1,
                "photon_encounter_ordinal": 0,
                "photon_bdt_score": 0.90,
                "photon_bdt_tight_threshold": 0.80,
                "photon_iso_r04": 1.0,
                "photon_iso_r04_threshold": 4.0,
            },
            {
                "photon_et": 24.0,
                "photon_eta": 0.2,
                "photon_encounter_ordinal": 1,
                "photon_bdt_score": 0.91,
                "photon_bdt_tight_threshold": 0.80,
                "photon_iso_r04": 2.0,
                "photon_iso_r04_threshold": 4.0,
            },
        ]
        self.assertEqual(program.choose_leader(records), 1)

    def test_noncanonical_photon_cuts_recompute_without_a_recorded_witness(self) -> None:
        selection = RecoilSelection(region="A", photon_et_min=20.0)
        program = compile_recoil_selection(selection)
        self.assertIsNone(program.leader_branch)
        self.assertIsNone(selection.recorded_leader_branch)
        with self.assertRaisesRegex(ValueError, "offline recomputation"):
            _ = selection.leader_branch
        records = [
            {
                "photon_et": 19.0,
                "photon_eta": 0.1,
                "photon_encounter_ordinal": 0,
                "photon_bdt_score": 0.90,
                "photon_bdt_tight_threshold": 0.80,
                "photon_iso_r04": 1.0,
                "photon_iso_r04_threshold": 4.0,
            },
            {
                "photon_et": 24.0,
                "photon_eta": 0.2,
                "photon_encounter_ordinal": 1,
                "photon_bdt_score": 0.91,
                "photon_bdt_tight_threshold": 0.80,
                "photon_iso_r04": 2.0,
                "photon_iso_r04_threshold": 4.0,
            },
        ]
        self.assertEqual(program.choose_leader(records), 1)


if __name__ == "__main__":
    unittest.main()
