"""Local mechanical fixtures; no production inputs or historical runtime imports."""

from pathlib import Path
import copy
import dataclasses
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import uproot
import yaml

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import features, train, calibration, registry, inputs, augment, frontend, diagnostics


def write_tree(f, name, columns):
    f.mktree(
        name,
        {
            k: ("string" if np.asarray(v).dtype.kind in "US" else np.asarray(v).dtype)
            for k, v in columns.items()
        },
    )
    f[name].extend(columns)


def base_file(path):
    # Non-monotonic identities exercise keyed joins, including 64-bit precision.
    eh = np.array([5, 5, 7], dtype=np.uint64)
    el = np.array([9, 9, 11], dtype=np.uint64)
    ph = np.array([100, 300, 200], dtype=np.uint64)
    pl = np.array([2**63 + 3, 4, 5], dtype=np.uint64)
    keys = dict(event_hi=eh, event_lo=el, photon_hi=ph, photon_lo=pl)
    with uproot.recreate(path) as f:
        f["metadata"] = (
            "contract=PhotonJetTrees\ncollision_system=1\ndata_kind=1\nsimulation_role=0\n"
        )
        f["completion"] = "completion_status=complete\nretained_events=2\n"
        write_tree(
            f,
            "Sources",
            dict(
                source_hi=np.array([2**63 + 1], dtype=np.uint64),
                source_lo=np.array([22], dtype=np.uint64),
                completed=np.array([1], np.int32),
            ),
        )
        write_tree(
            f,
            "Events",
            dict(
                source_hi=np.full(2, 2**63 + 1, np.uint64),
                source_lo=np.full(2, 22, np.uint64),
                event_hi=np.array([7, 5], np.uint64),
                event_lo=np.array([11, 9], np.uint64),
                centrality_percent=np.array([np.nan, np.nan]),
                centrality_valid=np.zeros(2, np.int32),
                reco_vertex_z=np.array([0.0, 1.0]),
                reco_vertex_valid=np.ones(2, np.int32),
                photon_count=np.array([1, 2], np.int32),
            ),
        )
        write_tree(
            f,
            "Photons",
            dict(
                keys,
                encounter_ordinal=np.arange(3, dtype=np.int64),
                energy=np.array([20.0, 25.0, 40.0]),
                et=np.array([20.0, 25.0, 40.0]),
                eta=np.array([0.1, 0.2, 0.1]),
                phi=np.zeros(3),
                kinematics_finite=np.ones(3, np.int32),
                producer_vertex_z=np.zeros(3),
            ),
        )
        shower = {k: np.full(3, 1.0, dtype=float) for k in features.SHOWER_SCALARS}
        shower["e33"] = np.full(3, 3.0)
        shower["e35"] = np.full(3, 5.0)
        order = [2, 0, 1]
        write_tree(
            f,
            "PhotonShowerViews",
            {
                k: np.asarray(v)[order]
                for k, v in dict(
                    keys, definition=np.array(["H70"] * 3), **shower
                ).items()
            },
        )
        write_tree(
            f,
            "Isolation",
            dict(
                keys,
                method=np.full(3, 3, np.int32),
                radius=np.full(3, 0.4),
                cone_sum=np.array([0.5, 4.0, 0.0]),
                valid=np.ones(3, np.int32),
                candidate_removed=np.ones(3, np.int32),
            ),
        )
    return path


def wp(system="pp"):
    return dict(
        schema="PhotonIDWorkingPointsV1",
        system=system,
        axis="et" if system == "pp" else "centrality",
        domain=[15, 35] if system == "pp" else [0, 80],
        fit_population="validation",
        curves={"WP70": [0.7, 0.0], "WP80": [0.6, 0.0], "WP90": [0.4, 0.0]},
        non_tight_lower=0.1,
        isolation=False,
    )


class ProtocolTests(unittest.TestCase):
    def test_group_union_and_order_independence(self):
        records = [
            ["hard1", "mb1"],
            ["hard2", "mb1"],
            ["hard2", "mb2"],
            ["hard3", "mb3"],
        ]
        groups = inputs.grouped_events(records)
        self.assertEqual(len(set(groups[:3])), 1)
        profile = registry.read(HERE / "training_profiles.yaml")["profiles"][
            "canonical_auau_v1"
        ]["split"]
        p = train.split_rows(groups, profile)
        perm = [3, 1, 2, 0]
        other = train.split_rows(
            inputs.grouped_events([records[i] for i in perm]), profile
        )
        np.testing.assert_array_equal(p[perm], other)
        self.assertEqual(len(set(p[:3])), 1)

    def test_group_partition_fractions_and_no_overlap(self):
        groups = np.repeat([f"event{i}" for i in range(20000)], 2)
        s = train.split_rows(groups, dict(rule="physical_event_hash_v1", seed=13))
        np.testing.assert_array_equal(s[::2], s[1::2])
        for name, target in [("train", 0.8), ("validation", 0.1), ("test", 0.1)]:
            self.assertAlmostEqual(np.mean(s == name), target, delta=0.012)

    def test_parity_split_exact_indices(self):
        y = np.tile([0, 1], 100)
        from sklearn.model_selection import train_test_split

        rest, test = train_test_split(
            np.arange(200), test_size=0.2, random_state=42, stratify=y
        )
        tr, val = train_test_split(
            rest, test_size=0.1 / (0.7 + 0.1), random_state=43, stratify=y[rest]
        )
        for a, b in zip(train.parity_split_indices(y), (tr, val, test)):
            np.testing.assert_array_equal(a, b)

    def test_flatten_per_source_deterministic(self):
        et = np.tile(np.linspace(6, 25, 501), 2)
        samples = np.repeat(["Jet8", "Jet12"], 501)
        y = np.zeros(1002, int)
        order = train.parity_row_order(
            et, y, samples, np.ones(1002, bool), ["Jet8", "Jet12"]
        )
        a = order[samples[order] == "Jet8"]
        b = order[samples[order] == "Jet12"] - 501
        np.testing.assert_array_equal(a, b)
        self.assertTrue(set(np.flatnonzero(et >= 15)).issubset(set(order)))
        self.assertEqual(len(set(order)), len(order))

    def test_weight_fit_does_not_see_heldout(self):
        rng = np.random.default_rng(4)
        et = rng.uniform(15, 35, 500)
        eta = rng.uniform(-0.7, 0.7, 500)
        y = np.tile([0, 1], 250)
        mask = np.arange(500) < 400
        config = dict(eta_range=[-0.7, 0.7], et_cap=800.0)
        a = train.fit_weights(et, eta, y, mask, config)
        et[400:] = 1e6
        eta[400:] = 1e6
        self.assertEqual(a, train.fit_weights(et, eta, y, mask, config))
        weights = train.apply_weights(et[:400], eta[:400], y[:400], a)
        self.assertTrue(np.all(np.isfinite(weights) & (weights > 0)))

    def test_weight_formula_matches_recovered_spline(self):
        from scipy.interpolate import UnivariateSpline

        rng = np.random.default_rng(11)
        et = rng.uniform(6, 35, 600)
        eta = rng.uniform(-0.7, 0.7, 600)
        y = np.tile([0, 1], 300)
        fitted = train.fit_weights(
            et, eta, y, np.ones(600, bool), dict(eta_range=[-0.7, 0.7], et_cap=800.0)
        )
        expected = np.ones(600)
        for c in (0, 1):
            m = y == c
            for x, bounds, cap in [
                (eta, (-0.7, 0.7), None),
                (et, (et[m].min(), et[m].max()), 800.0),
            ]:
                edges = np.linspace(*bounds, 21)
                h = np.histogram(x[m], bins=edges, density=True)[0] * 20
                w = 1 / np.clip(
                    UnivariateSpline((edges[:-1] + edges[1:]) / 2, h, s=0)(x[m]),
                    1e-3,
                    None,
                )
                if cap:
                    w = np.minimum(w, cap)
                expected[m] *= w / w.mean()
        np.testing.assert_allclose(
            train.apply_weights(et, eta, y, fitted), expected, rtol=1e-14
        )

    def test_profile_bindings_and_feature_order(self):
        specs = features.load_registry(HERE / "model_registry.yaml")
        profiles = registry.read(HERE / "training_profiles.yaml")["profiles"]
        for name, p in profiles.items():
            with self.assertRaisesRegex(ValueError, "unresolved"):
                train.check_profile(specs[name], p)
        self.assertEqual(
            [
                len(specs[n].features)
                for n in ("canonical_pp_v1", "canonical_auau_v1", "ppg12_equivalent_v1")
            ],
            [13, 14, 11],
        )


class CalibrationTests(unittest.TestCase):
    def test_exact_id_boundaries(self):
        s = np.array([np.nan, -1, 0.1, 0.10001, 0.599, 0.6, 0.7, 0.70001])
        np.testing.assert_array_equal(
            calibration.id_regions(s, 0.7, 0.6), [0, 1, 1, 2, 2, 3, 3, 4]
        )
        self.assertEqual(calibration.id_regions(0.8, 0.5, 0.6), 0)

    def test_exact_isolation_boundaries(self):
        np.testing.assert_array_equal(
            calibration.isolation_regions([np.nan, 0, 1, 1.5, 2, 3], 1, 2),
            [0, 1, 2, 2, 2, 3],
        )
        self.assertEqual(calibration.isolation_regions(1, 1, 1), 2)

    def test_weighted_quantile_ties(self):
        from numerics import weighted_threshold

        self.assertEqual(
            weighted_threshold(
                np.array([0.1, 0.1, 0.9]), np.array([1.0, 3.0, 1.0]), 0.7
            ),
            0.1,
        )

    def test_fit_only_validation(self):
        with self.assertRaisesRegex(ValueError, "validation"):
            calibration.fit_working_points([], [], [], [], [], "pp", partition="test")

    def test_fits_and_frozen_efficiency(self):
        et = np.repeat([16, 18, 20, 22, 24, 30], 200)
        score = np.tile(np.linspace(0.4, 0.99, 200), 6)
        w = np.ones(1200)
        cent = np.zeros(1200)
        p = calibration.fit_working_points(
            score, w, et, cent, w == 1, "pp", partition="validation"
        )
        calibration.validate_ordering(p)
        report = calibration.efficiency_report(score, w, et, cent, w == 1, p)
        self.assertTrue(
            all(abs(r["efficiency"] - int(r["wp"][2:]) / 100) < 0.015 for r in report)
        )
        before = registry.digest(p)
        calibration.efficiency_report(score * 0.9, w, et, cent, w == 1, p)
        self.assertEqual(before, registry.digest(p))

    def test_auau_isolation_grid(self):
        from itertools import product

        et = []
        cent = []
        v = []
        for c, p in product(np.arange(2.5, 80, 5), [16, 18, 20, 22, 24, 30]):
            et.extend([p] * 30)
            cent.extend([c] * 30)
            v.extend(np.linspace(-2, 10, 30) - 0.03 * c)
        n = len(v)
        p = calibration.fit_working_points(
            np.array(v),
            np.ones(n),
            np.array(et),
            np.array(cent),
            np.ones(n, bool),
            "auau",
            partition="validation",
            isolation=True,
        )
        self.assertEqual(len(p["tables"]["cells"]), 288)
        calibration.validate_ordering(p)

    def test_curve_crossing_refused(self):
        p = wp()
        p["curves"]["WP80"] = [0.8, 0]
        with self.assertRaises(ValueError):
            calibration.validate_ordering(p)

    def test_isolation_histogram_flow_refused(self):
        with self.assertRaisesRegex(ValueError, "overflow"):
            calibration.point(np.array([99.0, 100.0]), np.ones(2), 0.9, True)


class SidecarTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)
        self.base = base_file(self.d / "base.root")
        self.original = registry.file_hash(self.base)
        self.spec = features.load_registry(HERE / "model_registry.yaml")[
            "canonical_pp_v1"
        ]

    def bound_spec(self, name):
        artifact = self.d / (name + ".json")
        artifact.write_text("{}")
        receipt = self.d / (name + ".receipt.json")
        receipt.write_text("{}")
        return dataclasses.replace(
            self.spec,
            name=name,
            model_file=artifact,
            model_sha256=registry.file_hash(artifact),
            training_manifest=receipt,
        )

    def test_feature_join_and_identity_precision(self):
        m = features.features_for_file(self.base, self.spec)
        self.assertEqual(int(m.source[0, 0]), 2**63 + 1)
        self.assertEqual(m.keys[0][3], 2**63 + 3)
        np.testing.assert_allclose(
            m.values[:, 7], 1 / 3
        )  # canonical 13-feature ratio position
        np.testing.assert_array_equal(m.in_domain, [1, 1, 0])

    def test_multimodel_selection_root_roundtrip(self):
        specs = [self.bound_spec("model1"), self.bound_spec("model2")]
        iso = registry.read(HERE / "model_registry.yaml")["isolation_packages"][
            "ppg12_topocluster_r04_v1"
        ]
        recipe = dict(
            name="test_pp",
            model="model1",
            wp=wp(),
            isolation=iso,
            wp_sha256="a" * 64,
            isolation_sha256=registry.digest(iso),
        )
        source = dict(
            _path=self.base,
            sha256=self.original,
            source_id=[2**63 + 1, 22],
            lane="pp_data",
        )
        out = self.d / "sidecar.root"
        evaluators = {s.name: lambda x: np.full(len(x), 0.8) for s in specs}
        result = augment.augment_file(
            source,
            out,
            specs,
            evaluators,
            [recipe],
            {s.name: dict(version="1") for s in specs},
            {"name": "tiny_fixture"},
        )
        with uproot.open(out) as f:
            self.assertEqual(f["Models"].num_entries, 2)
            self.assertEqual(f["PhotonScores"].num_entries, 6)
            self.assertEqual(f["PhotonSelections"].num_entries, 3)
            np.testing.assert_array_equal(
                f["PhotonSelections"]["id_region"].array(library="np"), [4, 4, 0]
            )
            np.testing.assert_array_equal(
                f["PhotonSelections"]["isolation_region"].array(library="np"), [1, 3, 0]
            )
            self.assertEqual(
                int(f["PhotonScores"]["source_hi"].array(library="np")[0]), 2**63 + 1
            )
            self.assertEqual(f["PhotonScores"].classname, "TTree")
        self.assertEqual(self.original, registry.file_hash(self.base))
        self.assertEqual(result["sha256"], registry.file_hash(out))
        with self.assertRaises(FileExistsError):
            augment.augment_file(
                source,
                out,
                specs,
                evaluators,
                [recipe],
                {s.name: dict(version="1") for s in specs},
                {},
            )

    def test_raw_scores_need_no_recipe(self):
        s = self.bound_spec("raw")
        source = dict(
            _path=self.base,
            sha256=self.original,
            source_id=[2**63 + 1, 22],
            lane="pp_data",
        )
        augment.augment_file(
            source,
            self.d / "raw.root",
            [s],
            {"raw": lambda x: np.full(len(x), 0.5)},
            [],
            {"raw": dict(version="1")},
            {},
        )
        with uproot.open(self.d / "raw.root") as f:
            self.assertEqual(f["SelectionRecipes"].num_entries, 0)

    def test_invalid_feature_state_and_nonfinite_score(self):
        m = features.features_for_file(self.base, self.spec)
        m.complete[0] = False
        score, state = augment.score_matrix(m, lambda x: np.full(len(x), np.nan))
        np.testing.assert_array_equal(state, [3, 2, 4])
        self.assertTrue(np.all(np.isnan(score)))

    def test_event_source_identity_must_match_manifest(self):
        with uproot.open(self.base) as f:
            events = f["Events"].arrays(library="np")
        events["source_lo"][0] = 999
        with uproot.update(self.base) as f:
            write_tree(f, "Events", events)
        spec = self.bound_spec("raw")
        source = dict(
            _path=self.base,
            sha256=registry.file_hash(self.base),
            source_id=[2**63 + 1, 22],
            lane="pp_data",
        )
        with self.assertRaisesRegex(ValueError, "event source identity"):
            augment.augment_file(
                source,
                self.d / "invalid.root",
                [spec],
                {"raw": lambda x: np.full(len(x), 0.5)},
                [],
                {},
                {},
            )

    def test_wrong_source_and_mutated_base_rejected(self):
        s = self.bound_spec("raw")
        source = dict(
            _path=self.base, sha256=self.original, source_id=[1, 22], lane="pp_data"
        )
        with self.assertRaisesRegex(ValueError, "source mismatch"):
            augment.augment_file(source, self.d / "wrong.root", [s], {}, [], {}, {})
        source["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "base changed"):
            augment.augment_file(source, self.d / "wrong.root", [s], {}, [], {}, {})

    def test_model_wp_mismatch_rejected(self):
        spec = self.bound_spec("canonical_pp_v1")
        p = wp()
        p.update(
            model_sha256="b" * 64, feature_schema_sha256=features.feature_identity(spec)
        )
        path = self.d / "wp.json"
        registry.write_json(path, p)
        raw = registry.read(HERE / "model_registry.yaml")
        raw["selection_recipes"]["test"] = dict(
            model=spec.name,
            id_package="wp",
            isolation_package="ppg12_topocluster_r04_v1",
        )
        raw["working_points"]["wp"] = dict(
            path=str(path), sha256=registry.file_hash(path)
        )
        with self.assertRaisesRegex(ValueError, "different model"):
            registry.resolve_selection(raw, "test", {spec.name: spec}, self.d)
        p["model_sha256"] = spec.model_sha256
        registry.write_json(path, p)
        raw["working_points"]["wp"]["sha256"] = registry.file_hash(path)
        self.assertEqual(
            registry.resolve_selection(raw, "test", {spec.name: spec}, self.d)[
                "wp_sha256"
            ],
            registry.file_hash(path),
        )

    def test_prompt_labels_are_independent_of_analysis_isolation(self):
        with uproot.open(self.base) as f:
            photons = f["Photons"].arrays(library="np")
            events = f["Events"].arrays(library="np")
        photons.update(
            dominant_truth_state=np.array([2, 2, 2], np.int32),
            dominant_truth_pid=np.array([22, 22, 211], np.int32),
            dominant_truth_track_id=np.array([1, 2, 3], np.int32),
            dominant_truth_embedding_id=np.array([0, 0, 0], np.int32),
        )
        events["truth_denominator_complete"] = np.ones(2, np.int32)
        with uproot.update(self.base) as f:
            write_tree(f, "Photons", photons)
            write_tree(f, "Events", events)
            write_tree(
                f,
                "TruthPhotons",
                dict(
                    event_hi=np.array([5, 5], np.uint64),
                    event_lo=np.array([9, 9], np.uint64),
                    track_id=np.array([1, 2], np.int32),
                    embedding_id=np.array([0, 0], np.int32),
                    prompt_class=np.array([1, 3], np.int32),
                    generator_association_valid=np.ones(2, np.int32),
                    analysis_signal=np.array([0, 1], np.int32),
                ),
            )
        m = features.features_for_file(self.base, self.spec)
        label, analysis = inputs.dominant_prompt_labels(self.base, m)
        np.testing.assert_array_equal(label, [1, 0, 0])
        np.testing.assert_array_equal(analysis, [False, True, False])

    def test_receipt_model_binding_is_checked(self):
        spec = self.bound_spec("canonical_pp_v1")
        raw = registry.read(HERE / "model_registry.yaml")
        r = raw["models"][spec.name]
        r.update(
            status="candidate",
            binding_evidence="synthetic test",
            model_file=str(spec.model_file),
            model_sha256=spec.model_sha256,
            feature_schema_sha256=features.feature_identity(spec),
            training_manifest=str(spec.training_manifest),
            training_receipt_sha256=registry.file_hash(spec.training_manifest),
        )
        path = self.d / "registry.yaml"
        path.write_text(yaml.safe_dump(raw))
        with self.assertRaisesRegex(ValueError, "receipt belongs"):
            features.load_registry(path)
        registry.write_json(
            spec.training_manifest,
            dict(
                model_sha256=spec.model_sha256,
                feature_schema_sha256=features.feature_identity(spec),
            ),
        )
        r["training_receipt_sha256"] = registry.file_hash(spec.training_manifest)
        path.write_text(yaml.safe_dump(raw))
        self.assertTrue(features.load_registry(path)[spec.name].bound)

    def test_xgboost_serialization_evaluator(self):
        from xgboost import XGBClassifier
        import pandas as pd

        rng = np.random.default_rng(7)
        x = rng.normal(size=(100, 13)).astype(np.float32)
        y = (x[:, 0] > 0).astype(int)
        m = XGBClassifier(n_estimators=3, max_depth=2, n_jobs=1)
        m.fit(pd.DataFrame(x, columns=self.spec.features), y)
        path = self.d / "tiny.json"
        m.save_model(path)
        spec = dataclasses.replace(
            self.spec, model_file=path, model_sha256=registry.file_hash(path)
        )
        np.testing.assert_array_equal(
            train.predict(m, x, spec.features), augment.make_evaluator(spec)(x)
        )


class FrontendTests(unittest.TestCase):
    def test_composition(self):
        self.assertEqual(registry.expand("runAll"), list(registry.LANES))
        self.assertEqual(registry.expand("runAllData"), ["runPP", "runAuAu"])

    def test_help_from_arbitrary_directory(self):
        for name in ("train", "augment"):
            for arg in ([], ["--help"]):
                r = subprocess.run(
                    [str(HERE / f"{name}_photon_id.sh"), *arg],
                    cwd="/tmp",
                    text=True,
                    capture_output=True,
                )
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertIn("sPHENIX Photon ID", r.stdout)

    def test_unknown_command_and_fail_closed(self):
        env = dict(os.environ, PHOTONID_PYTHON=sys.executable)
        for script, command in [
            ("train", "trainAll"),
            ("augment", "runAll"),
            ("train", "bogus"),
            ("augment", "bogus"),
        ]:
            r = subprocess.run(
                [str(HERE / f"{script}_photon_id.sh"), command, "--dry-run"],
                env=env,
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("ERROR:", r.stderr)

    def test_incompatible_model_rejected(self):
        specs = features.load_registry(HERE / "model_registry.yaml")
        args = frontend.parser("augment").parse_args(
            ["runAuAu", "--models", "canonical_pp_v1"]
        )
        with self.assertRaisesRegex(ValueError, "incompatible"):
            frontend.model_choices(
                "auau", args, registry.read(HERE / "model_registry.yaml"), specs
            )

    def test_missing_campaign_and_duplicate_sources(self):
        with self.assertRaisesRegex(ValueError, "unbound"):
            registry.campaign(None, None)
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            f = d / "base"
            f.write_text("fixture")
            s = dict(
                path="base",
                sha256=registry.file_hash(f),
                lane="pp_data",
                source_id=[1, 2],
            )
            m = d / "manifest.json"
            registry.write_json(
                m,
                dict(
                    schema="PhotonJetTreeCampaignV1",
                    status="complete",
                    campaign="tiny",
                    sources=[s, s],
                ),
            )
            with self.assertRaisesRegex(ValueError, "duplicate"):
                registry.campaign(m, registry.file_hash(m))


class PackageTests(unittest.TestCase):
    def test_tiny_real_training_freeze_then_render(self):
        # Intentional tiny synthetic training, not an approved production profile.
        spec = features.load_registry(HERE / "model_registry.yaml")["canonical_pp_v1"]
        profile = copy.deepcopy(
            registry.read(HERE / "training_profiles.yaml")["profiles"][spec.name]
        )
        profile["xgboost"].update(n_estimators=25, max_depth=3, n_jobs=1)
        rng = np.random.default_rng(12)
        n = 3000
        et = rng.uniform(15, 35, n)
        eta = rng.uniform(-0.7, 0.7, n)
        X = rng.normal(size=(n, 13)).astype(np.float32)
        X[:, 0] = et
        y = (rng.uniform(size=n) < 1 / (1 + np.exp(-X[:, 1] * 0.8))).astype(int)
        data = dict(
            X=X,
            y=y,
            et=et,
            eta=eta,
            centrality=np.zeros(n),
            ready=np.ones(n, bool),
            analysis=y == 1,
            isolation=rng.normal(size=n),
            groups=np.array([f"e{i//2}" for i in range(n)]),
            keys=np.column_stack(
                (np.ones((n, 5), np.uint64), np.arange(n, dtype=np.uint64))
            ),
            sample=np.where(y == 1, "Photon20", "Jet20"),
        )
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "package"
            original = train.predict
            calls = []

            def spy(model, values, names):
                calls.append((out / "FROZEN_CHOICES.json").exists())
                return original(model, values, names)

            with patch.object(train, "predict", side_effect=spy):
                receipt = train.train_package(
                    spec, profile, data, out, dict(fixture=True)
                )
            self.assertEqual(calls, [False, False, True])
            self.assertEqual(receipt["mechanical_status"], "PASS")
            self.assertEqual(receipt["scientific_acceptance"], "NOT_REVIEWED")
            self.assertFalse((out / "INCOMPLETE.json").exists())
            for population in ("validation", "test", "all"):
                directory = out / "diagnostics" / population
                self.assertEqual(len(list(directory.glob("*.png"))), 4)
                report = registry.read(directory / "numerical.json")
                self.assertEqual(report["wp"]["fit_population"], "validation")
            self.assertEqual(
                receipt["metrics"]["all"]["label"], "full-sample / not held-out"
            )
            self.assertIn("split/weighting.json", receipt["files"])
            for name, h in receipt["files"].items():
                self.assertEqual(registry.file_hash(out / name), h)


if __name__ == "__main__":
    unittest.main()
