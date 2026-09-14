"""Synthetic transport/boundary tests; never a real-DST acceptance receipt."""
from pathlib import Path
import copy
import json
import math
import os
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import awkward as ak
import numpy as np
import uproot

import centrality_replay as replay
import build_photonjet_collaboration_tree as builder
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "TreeProduction" / "producer"))
import prepare_centrality_source_delta as overlay


NODES = dict(pmt="MbdPmtContainer", mbd_out="MbdOut", minimum_bias="MinimumBiasInfo", centrality="CentralityInfo")


def calibration():
    return dict(divisions=np.arange(100, 0, -1, dtype="float32"), run_scale=1.0,
                vertex_ranges=[[-60., 0., 1.], [0., 60., 1.25]])


def events():
    n = 5
    x = {k: np.full(n, -1 if v == "int32" else np.nan, dtype=v) for k, v in replay.SCALARS.items()}
    for k in ("centrality_replay_version", "centrality_pmt_available", "centrality_inputs_valid",
              "centrality_mbd_z_valid", "centrality_native_valid"):
        x[k][:] = 1
    x.update({k: np.arange(n, dtype="uint64") + j + 1 for j, k in enumerate(replay.KEYS)})
    x.update(run=np.full(n, 68390, dtype="int32"), centrality=np.full(n, 50.),
             centrality_mb_decision=np.asarray([1, 1, 0, 1, 1], dtype="int32"),
             centrality_mbd_z=np.asarray([0., 60., 0., -60., np.nan]),
             centrality_selected_charge=np.full(n, 64.))
    x["centrality_mbd_z_valid"][-1] = 0
    x["centrality_inputs_valid"][-1] = 0
    for k, row in dict(centrality_pmt_id=list(range(128)), centrality_pmt_charge=[.5] * 128,
                       centrality_pmt_time=[25.] * 128, centrality_pmt_valid=[1] * 128,
                       centrality_pmt_selected=[1] * 128).items():
        x[k] = ak.Array([row for _ in range(n)])
    return x


def root_file(path, tree, schema, arrays):
    with uproot.recreate(path) as f:
        f.mktree(tree, schema).extend(arrays)


class ReplayTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.home = Path(self.tmp.name)

    def manifest(self):
        """Build genuine CDBTTree-shaped ROOT fixtures, explicitly synthetic."""
        payloads = {}
        for domain, fields in {
            "Centrality": {"IID": np.arange(100, dtype="int32"), "Fcentralitydiv": calibration()["divisions"]},
            "CentralityScale": {"IID": np.array([0], dtype="int32"), "Dcentralityscale": np.array([1.])},
            "CentralityVertexScale": {"IID": np.array([0, 1], dtype="int32"),
                "Invertexbins": np.array([2, -2147483648], dtype="int32"),
                "Dscale": np.array([1., 1.25]), "Dlow_vertex": np.array([-60., 0.]), "Dhigh_vertex": np.array([0., 60.])},
        }.items():
            path = self.home / (domain + ".root")
            root_file(path, "Multiple", {k: v.dtype for k, v in fields.items()}, fields)
            payloads[domain] = dict(path=str(path), sha256=replay.digest(path))
        software = {"release": "SYNTHETIC_FIXTURE_ONLY"}
        hashes = {}
        for name in replay.SUPPORTED_NATIVE_SHA256:
            path = self.home / name
            path.write_text("// Synthetic software identity fixture, no production authority\n")
            hashes[name] = replay.digest(path)
            software[name] = dict(path=str(path), sha256=hashes[name])
        self.addCleanup(patch.stopall)
        patch.object(replay, "SUPPORTED_NATIVE_SHA256", hashes).start()
        self.manifest_data = dict(schema="AuAuCentralityCalibrationManifestV1", algorithm=replay.ALGORITHM,
            software=software, records=[dict(run=68390, global_tag="SYNTHETIC_FIXTURE_ONLY", nodes=NODES, payloads=payloads)])
        p = self.home / "manifest.json"
        p.write_text(json.dumps(self.manifest_data))
        return p

    def inputs(self, data=None):
        self.data = events() if data is None else data
        source = self.home / "source.root"
        schema = {**{k: v.dtype for k, v in self.data.items() if k not in replay.ARRAYS}, **replay.ARRAYS}
        root_file(source, "eventTree", schema, self.data)
        self.binding_data = dict(schema="AuAuCentralityInputBindingV1", source_sha256=replay.digest(source),
            system="AuAu", sample_kind="data", nodes=NODES, calibration_run_branch="run",
            input_measurement_definition=replay.ALGORITHM)
        p = self.home / "binding.json"
        p.write_text(json.dumps(self.binding_data))
        return source, p

    def test_cdb_payloads_and_exact_boundaries(self):
        _, c = replay.load_calibrations(self.manifest())
        self.assertEqual(c[68390]["vertex_ranges"], calibration()["vertex_ranges"])
        r = replay.calculate([.5] * 128, [25] * 128, [1] * 128, 0, 1, c[68390])
        self.assertEqual((r["scaled_charge"], r["bin"]), (64., 38))  # division == charge is skipped
        self.assertEqual(replay.calculate([.5] * 128, [-25] * 128, [1] * 128, 60, 1, c[68390])["scaled_charge"], 80.)
        self.assertEqual(replay.calculate([.5] * 128, [25] * 128, [1] * 128, -60, 1, c[68390])["state"], 3)

    def test_roundtrip_keeps_all_rows_and_originals(self):
        m = self.manifest(); src, b = self.inputs()
        out = self.home / "friend"
        r = replay.augment(src, m, b, out, step_size=2)
        self.assertEqual(r["states"], dict(VALID=2, NON_MB=1, VERTEX_OUTSIDE_PAYLOAD=1, MISSING_OR_INVALID_INPUT=1))
        self.assertFalse(r["production_accepted"])
        with uproot.open(out / "centrality_friend.root") as f:
            values = f["centralityFriend"].arrays(library="np")
        np.testing.assert_equal(values["source_entry"], np.arange(5))
        for key in replay.KEYS:
            np.testing.assert_equal(values[key], self.data[key])
        np.testing.assert_equal(values["centrality_original"], self.data["centrality"])
        self.assertEqual(values["centrality_bin_new"].tolist(), [38, 22, 0, 0, 0])
        with self.assertRaisesRegex(ValueError, "already exists"):
            replay.augment(src, m, b, out)

    def test_duplicate_occurrence_rejected_across_chunks(self):
        m = self.manifest(); x = events()
        for key in replay.KEYS:
            x[key][-1] = x[key][0]
        src, b = self.inputs(x)
        with self.assertRaisesRegex(Exception, "UNIQUE constraint"):
            replay.augment(src, m, b, self.home / "bad", step_size=1)
        self.assertEqual(json.loads((self.home / "bad/RECEIPT.json").read_text())["status"], "FAIL")

    def test_normalized_multi_source_keys_and_tampered_source_index(self):
        m = self.manifest()
        base = events()
        x = {k: (ak.concatenate([v,v]) if k in replay.ARRAYS else np.concatenate([v,v]))
             for k,v in base.items()}
        x["source_file_index"] = np.repeat(np.array([0,1], dtype="int32"), 5)
        src, b = self.inputs(x)
        # The normalized layout uses events, and compact references repeat per source.
        schema = {**{k:v.dtype for k,v in x.items() if k not in replay.ARRAYS}, **replay.ARRAYS}
        root_file(src, "events", schema, x)
        self.binding_data["source_sha256"] = replay.digest(src)
        b.write_text(json.dumps(self.binding_data))
        out = self.home / "normalized_friend"
        receipt = replay.augment(src, m, b, out, step_size=3)
        self.assertEqual(receipt["key_fields"], [*replay.KEYS,"source_file_index"])
        self.assertEqual(receipt["source_tree"], "events")
        self.assertEqual(replay.verify(src, out, step_size=4)["rows"], 10)
        root_path = out / "centrality_friend.root"
        with uproot.open(root_path) as f:
            values = f["centralityFriend"].arrays(library="np")
        values["source_file_index"][-1] = 2  # no duplicate, but incorrect join target
        root_file(root_path,"centralityFriend",{k:v.dtype for k,v in values.items()},values)
        receipt["friend_sha256"] = replay.digest(root_path)
        (out / "RECEIPT.json").write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,"has no friend row"):
            replay.verify(src,out)

    def test_wrong_run_or_embedded_identity_rejected(self):
        m = self.manifest(); x = events(); x["run"][-1] = 68144
        src, b = self.inputs(x)
        with self.assertRaisesRegex(ValueError, "exact calibration binding"):
            replay.augment(src, m, b, self.home / "bad")
        self.binding_data["sample_kind"] = "embedded"
        b.write_text(json.dumps(self.binding_data))
        with self.assertRaisesRegex(ValueError, "background-event identity"):
            replay.augment(src, m, b, self.home / "bad2")

    def test_payload_and_source_drift_rejected(self):
        m = self.manifest(); src, b = self.inputs()
        with src.open("ab") as f:
            f.write(b"changed")
        with self.assertRaisesRegex(ValueError, "binding/source hash"):
            replay.augment(src, m, b, self.home / "bad")
        spec = self.manifest_data["records"][0]["payloads"]["CentralityScale"]
        Path(spec["path"]).write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "payload hash"):
            replay.load_calibrations(m)

    def test_invalid_and_no_division(self):
        c = calibration()
        for q, t, mask, z, mb, expected in [
            ([np.nan] * 128, [0.] * 128, [1] * 128, 0, 1, 2),
            ([1.] * 128, [np.nan] * 128, [1] * 128, 0, 1, 2),
            ([1.] * 128, [0.] * 128, [0] * 128, 0, 1, 2),
            ([.4999] * 128, [0.] * 128, [1] * 128, 0, 1, 4),
            ([1.] * 128, [25.001] * 128, [1] * 128, 0, 1, 4),
            ([1.] * 128, [0.] * 128, [1] * 128, 0, -1, 2),
        ]:
            self.assertEqual(replay.calculate(q, t, mask, z, mb, c)["state"], expected)

    def test_channel_and_mask_corruption(self):
        for field, value, message in (("centrality_pmt_id", 99, "channel identity"),
                                      ("centrality_pmt_selected", 0, "PMT selection")):
            with self.subTest(field=field):
                x = events(); a = ak.to_list(x[field]); a[0][0] = value; x[field] = ak.Array(a)
                m = self.manifest(); src, b = self.inputs(x)
                with self.assertRaisesRegex(ValueError, message):
                    replay.augment(src, m, b, self.home / field)

    def test_legacy_export_and_partial_capture(self):
        e = {"event_id_hi": np.array([1], dtype="uint64")}
        s, a = replay.capture_columns(e)
        self.assertEqual(s["centrality_replay_version"][0], 0)
        self.assertEqual(ak.to_list(a["centrality_pmt_time"]), [[]])
        e["centrality_replay_version"] = np.array([1])
        with self.assertRaisesRegex(ValueError, "partial"):
            replay.capture_columns(e)

    def test_export_preserves_exact_replay_arrays(self):
        x = events(); x["vertex_z"] = x["centrality_mbd_z"]
        s = builder._interface_columns(x, None, "auau")
        for k in replay.SCALARS:
            np.testing.assert_equal(s[k], x[k])
        _, a = replay.capture_columns(x)
        self.assertEqual(ak.to_list(a["centrality_pmt_time"]), ak.to_list(x["centrality_pmt_time"]))

    def test_complete_exporter_root_roundtrip_with_replay_capture(self):
        import test_collaborator_interface as interface
        import validate_photonjet_collaboration_tree as validator
        e = interface.event_fixture()
        n = len(e["event_id_hi"])
        sample = events()
        for k in replay.SCALARS:
            e[k] = np.repeat(sample[k][:1], n)
        for k in replay.ARRAYS:
            e[k] = ak.Array([ak.to_list(sample[k][0])] * n)
        source, output = self.home / "internal.root", self.home / "collaborator.root"
        with patch.object(interface, "event_fixture", return_value=e):
            interface.write_source(source, system="auau")
        builder.build([source], output, "auau", require_complete_interface=True, require_centrality_replay=True)
        with uproot.open(output) as root:
            validator.validate_interface([source], root, "auau", require_complete=True)
            values = root["eventTree"].arrays(list(replay.ARRAYS), library="ak")
            self.assertEqual(ak.to_list(values["centrality_pmt_time"]), [[25.] * 128] * n)
            np.testing.assert_equal(root["events"]["centrality_mbd_z"].array(library="np"), e["centrality_mbd_z"])

    def test_required_capture_rejects_charge_only_source(self):
        import test_collaborator_interface as interface
        source = self.home / "charge_only.root"
        interface.write_source(source, system="auau")
        with self.assertRaisesRegex(ValueError, "complete centrality replay capture"):
            builder.build([source], self.home / "export.root", "auau", require_centrality_replay=True)

    def test_capture_header_missing_timing_boundaries_and_reset(self):
        declarations = "\n".join(f"{'int' if t=='int32' else 'double'} {k};" for k, t in replay.SCALARS.items())
        declarations += "\n" + "\n".join(f"std::vector<{'int' if t.endswith('int32') else 'double'}> {k};" for k, t in replay.ARRAYS.items())
        code = '#include "RJCentralityReplayV1.h"\n#include <vector>\n#include <cassert>\n'
        code += 'struct Row {' + declarations + '};\n'
        code += '''
struct Hit {int id;float q=.5f,t=25;int get_pmt()const{return id;}float get_q()const{return q;}float get_time()const{return t;}};
struct Pmts {std::vector<Hit> hits;int missing=-1;int get_npmt(){return hits.size();}Hit* get_pmt(int i){return i==missing?nullptr:&hits[i];}};
struct Out {float z=60;float get_zvtx(){return z;}int get_evt(){return 123;}int get_clock(){return 321;}int get_femclock(){return 88;}};
struct MB {bool accepted=true;bool isAuAuMinimumBias(){return accepted;}};
int main(){Row r;Pmts p;Out o;MB m;for(int i=0;i<128;++i)p.hits.push_back({i});
RJCentralityReplayV1::capture(r,&p,&o,&m);
assert(r.centrality_inputs_valid==1 && r.centrality_selected_charge==64 && r.centrality_mbd_z==60);
assert(r.centrality_mbd_event==123 && r.centrality_mbd_clock==321 && r.centrality_pmt_time[0]==25);
p.hits[0].t=-25;p.hits[1].q=.499f;p.hits[2].t=25.001f;
RJCentralityReplayV1::capture(r,&p,&o,&m);assert(r.centrality_selected_charge==63);
assert(r.centrality_pmt_selected[0]==1 && r.centrality_pmt_selected[1]==0 && r.centrality_pmt_selected[2]==0);
p.hits[3].t=std::numeric_limits<float>::quiet_NaN();
RJCentralityReplayV1::capture(r,&p,&o,&m);assert(r.centrality_inputs_valid==0 && std::isnan(r.centrality_selected_charge));
p.hits[3].t=0;p.missing=4;RJCentralityReplayV1::capture(r,&p,&o,&m);assert(r.centrality_pmt_valid[4]==0);
p.missing=-1;p.hits[4].id=9;RJCentralityReplayV1::capture(r,&p,&o,&m);assert(r.centrality_pmt_valid[4]==0);
RJCentralityReplayV1::capture(r,(Pmts*)nullptr,(Out*)nullptr,(MB*)nullptr);
assert(r.centrality_pmt_charge.empty() && r.centrality_inputs_valid==0 && r.centrality_mb_decision==-1 && std::isnan(r.centrality_mbd_z));
p.hits.pop_back();RJCentralityReplayV1::capture(r,&p,&o,&m);assert(r.centrality_pmt_available==-1 && r.centrality_pmt_charge.empty());
}
'''
        cpp = self.home / "capture.cc"; cpp.write_text(code); binary = self.home / "capture"
        subprocess.run(["clang++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
            "-I", str(Path(__file__).resolve().parents[1] / "TreeProduction" / "producer"), str(cpp), "-o", str(binary)], check=True, capture_output=True)
        subprocess.run([str(binary)], check=True)

    def test_overlay_exact_hashes_and_explicit_embedding_nodes(self):
        root = self.home / "base"; (root / "src").mkdir(parents=True); (root / "src_AuAu").mkdir()
        header = root / "src/RJReplayFoundationV1.h"
        header.write_text('  double mbd_total_charge=std::numeric_limits<double>::quiet_NaN();\n'
                          'scalar(m_tEvent,"mbd_total_charge",&m_event.mbd_total_charge,"D");\n')
        emission = root / "src_AuAu/RecoilJets_AuAu.cc"
        emission.write_text('RJReplayRuntimeV1::EventBundle bundle;\n// retain original truth and weight capture\n')
        binding = dict(schema="CentralityCaptureSourceBindingV1", nodes=NODES, sample_kind="data",
                       source_root=str(root), emission_source=str(emission.relative_to(root)),
                       source_sha256={str(p.relative_to(root)): replay.digest(p) for p in (header, emission)})
        path = self.home / "overlay.json"; path.write_text(json.dumps(binding))
        result = overlay.prepare(path, self.home / "overlay")
        self.assertFalse(result["production_allowed"])
        changed = (self.home / "overlay/src_AuAu/RecoilJets_AuAu.cc").read_text()
        self.assertIn('getClass<MbdOut>(topNode,"MbdOut")', changed)
        self.assertIn("retain original truth and weight capture", changed)
        binding["sample_kind"] = "embedded"; path.write_text(json.dumps(binding))
        with self.assertRaisesRegex(ValueError, "background node"):
            overlay.prepare(path, self.home / "bad")
        binding["sample_kind"] = "data"; path.write_text(json.dumps(binding)); header.write_text("drift")
        with self.assertRaisesRegex(ValueError, "source drift"):
            overlay.prepare(path, self.home / "bad2")


def parity_events():
    """A source whose native witness is exactly what `calibration()` produces.

    Only events 0 and 1 have a valid native centrality; the non-MB, outside
    payload and missing-vertex rows must stay explicitly witnessless.
    """
    x = events()
    x["centrality_native_bin"] = np.asarray([38, 22, -1, -1, -1], dtype="int32")
    x["centrality_native_centile"] = np.asarray(
        [float(np.float32(np.float32(b) / np.float32(100))) if b > 0 else np.nan
         for b in x["centrality_native_bin"]])
    x["centrality_native_valid"] = np.asarray([1, 1, 0, 0, 0], dtype="int32")
    return x


class ParityAndConsumptionTests(ReplayTests):
    """Native-versus-offline parity and the friend join a consumer must trust."""

    def target_manifest(self, divisions):
        """A second, different calibration: the one we would migrate to."""
        data = copy.deepcopy(self.manifest_data)
        payloads = data["records"][0]["payloads"]
        path = self.home / "TargetCentrality.root"
        root_file(path, "Multiple", {"IID": np.dtype("int32"), "Fcentralitydiv": np.dtype("float32")},
                  {"IID": np.arange(100, dtype="int32"),
                   "Fcentralitydiv": np.asarray(divisions, dtype="float32")})
        payloads["Centrality"] = dict(path=str(path), sha256=replay.digest(path))
        target = self.home / "target_manifest.json"
        target.write_text(json.dumps(data))
        return target

    def test_native_parity_passes_against_the_as_produced_calibration(self):
        produced = self.manifest()
        src, b = self.inputs(parity_events())
        # Raise every division by half so the same charge lands in a new bin.
        target = self.target_manifest(np.arange(100, 0, -1) * 1.5)
        out = self.home / "friend"
        receipt = replay.augment(src, target, b, out, parity_manifest=produced)
        self.assertEqual(receipt["native_parity"], "PASS")
        self.assertEqual(receipt["parity_counts"], {"checked": 2, "agree": 2, "no_native_witness": 3})
        with uproot.open(out / "centrality_friend.root") as f:
            new_bins = f["centralityFriend"]["centrality_bin_new"].array(library="np")
        self.assertEqual(new_bins.tolist(), [59, 48, 0, 0, 0])  # native was 38/22

    def test_native_parity_rejects_a_disagreeing_witness(self):
        produced = self.manifest()
        x = parity_events()
        x["centrality_native_bin"][0] = 39  # one bin off
        src, b = self.inputs(x)
        target = self.target_manifest(np.arange(100, 0, -1) * 1.5)
        with self.assertRaisesRegex(ValueError, "parity failed"):
            replay.augment(src, target, b, self.home / "bad", parity_manifest=produced)
        receipt = json.loads((self.home / "bad/RECEIPT.json").read_text())
        self.assertEqual((receipt["status"], receipt["native_parity"]), ("FAIL", "FAIL"))

    def test_parity_manifest_must_be_a_different_calibration(self):
        produced = self.manifest()
        src, b = self.inputs(parity_events())
        with self.assertRaisesRegex(ValueError, "as-produced calibration"):
            replay.augment(src, produced, b, self.home / "same", parity_manifest=produced)

    def test_absent_parity_manifest_is_recorded_as_unchecked(self):
        m = self.manifest(); src, b = self.inputs(parity_events())
        out = self.home / "unchecked"
        receipt = replay.augment(src, m, b, out)
        self.assertEqual(receipt["native_parity"], "NOT_CHECKED")
        report = replay.verify(src, out)
        self.assertFalse(report["safe_to_consume_as_native_reproduction"])

    def test_verify_proves_the_join_and_reports_boundary_migration(self):
        produced = self.manifest()
        src, b = self.inputs(parity_events())
        # Raise every division so both events fall outside the 0-80% window.
        target = self.target_manifest(np.arange(100, 0, -1) + 1000.0)
        out = self.home / "friend"
        replay.augment(src, target, b, out, parity_manifest=produced, step_size=2)
        report = replay.verify(src, out, step_size=3)  # a step that splits differently
        self.assertEqual(report["rows"], 5)
        self.assertTrue(report["safe_to_consume_as_native_reproduction"])
        self.assertEqual(report["states"]["NO_DIVISION_MATCH"], 2)
        self.assertEqual(report["boundary"], {"left_analysis_window": 2})
        self.assertEqual(report["bin_changed"], 2)
        self.assertEqual(json.loads((out / "VERIFICATION.json").read_text())["rows"], 5)

    def test_verify_rejects_an_entry_aligned_friend_with_rotated_identities(self):
        m = self.manifest(); src, b = self.inputs(parity_events())
        out = self.home / "friend"
        receipt = replay.augment(src, m, b, out)
        root_path = out / "centrality_friend.root"
        with uproot.open(root_path) as f:
            values = f["centralityFriend"].arrays(library="np")
        for key in replay.KEYS:  # same row count, same order, wrong identities
            values[key] = np.roll(values[key], 1)
        schema = {k: v.dtype for k, v in values.items()}
        root_file(root_path, "centralityFriend", schema, values)
        receipt["friend_sha256"] = replay.digest(root_path)
        (out / "RECEIPT.json").write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "claims source entry"):
            replay.verify(src, out)

    def test_verify_rejects_a_mutated_preserved_original(self):
        m = self.manifest(); src, b = self.inputs(parity_events())
        out = self.home / "friend"
        receipt = replay.augment(src, m, b, out)
        root_path = out / "centrality_friend.root"
        with uproot.open(root_path) as f:
            values = f["centralityFriend"].arrays(library="np")
        values["centrality_original"][0] += 1.0
        root_file(root_path, "centralityFriend", {k: v.dtype for k, v in values.items()}, values)
        receipt["friend_sha256"] = replay.digest(root_path)
        (out / "RECEIPT.json").write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "preserved original differs"):
            replay.verify(src, out)

    def test_uncovered_vertex_domain_is_named_not_silently_accepted(self):
        """Native writes a bin from a zero estimator outside the vertex payload.

        `getVertexScale()` returns 0 when no interval matches, so FillVars sums
        zero and FillCentralityInfo still assigns a bin whenever any division is
        below zero.  The replay refuses to reproduce that; the operator has to
        acknowledge the divergence explicitly.
        """
        produced = self.manifest()
        x = parity_events()
        # Event 3 has z=-60, which no interval covers.  Give it the witness the
        # native pass would have written from a zero estimator.
        x["centrality_native_bin"][3] = 100
        x["centrality_native_centile"][3] = 1.0
        x["centrality_native_valid"][3] = 1
        src, b = self.inputs(x)
        target = self.target_manifest(np.arange(100, 0, -1) * 1.5)
        with self.assertRaisesRegex(ValueError, "uncovered vertex domain"):
            replay.augment(src, target, b, self.home / "strict", parity_manifest=produced)
        out = self.home / "acknowledged"
        receipt = replay.augment(src, target, b, out, parity_manifest=produced,
                                 allow_vertex_domain_divergence=True)
        self.assertEqual(receipt["native_parity"], "PASS_WITH_VERTEX_DOMAIN_DIVERGENCE")
        self.assertEqual(receipt["parity_counts"]["vertex_domain_divergence"], 1)
        report = replay.verify(src, out)
        self.assertFalse(report["safe_to_consume_as_native_reproduction"])
        self.assertIn("zero estimator", report["divergence"])

    def test_verify_rejects_a_friend_paired_with_another_source(self):
        m = self.manifest(); src, b = self.inputs(parity_events())
        out = self.home / "friend"
        replay.augment(src, m, b, out)
        with src.open("ab") as handle:
            handle.write(b"changed")
        with self.assertRaisesRegex(ValueError, "source changed"):
            replay.verify(src, out)


class NativeParityTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get("CENTRALITY_NATIVE_ROOT"), "set CENTRALITY_NATIVE_ROOT to pinned coresoftware source")
    def test_native_methods_bitwise(self):
        """Compile the actual three method bodies with minimal data-node fixtures.

        This is local source-algorithm parity, not a Fun4All/runtime/DST test.
        """
        source_root = Path(os.environ["CENTRALITY_NATIVE_ROOT"])
        source = source_root / "CentralityReco.cc"
        for name, expected in replay.SUPPORTED_NATIVE_SHA256.items():
            self.assertEqual(replay.digest(source_root / name), expected)
        text = source.read_text()
        def method(signature):
            start = text.index(signature); b = text.index("{", start); depth = 1; e = b + 1
            while depth:
                depth += (text[e] == "{") - (text[e] == "}"); e += 1
            return text[start:e]
        code = '''#include <array>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdint>
namespace Fun4AllReturnCodes {const int EVENT_OK=0;}
struct MbdPmtHit {float q,t; float get_q(){return q;} float get_time(){return t;}};
struct MbdPmtContainer {std::array<MbdPmtHit,128> pmts; MbdPmtHit* get_pmt(int i){return &pmts[i];}};
struct MbdOut {float z; float get_zvtx(){return z;}};
struct CentralityInfo {enum PROP{mbd_NS}; float c;int b;
 void set_centile(PROP,float v){c=v;} void set_centrality_bin(PROP,int v){b=v;}};
struct CentralityReco {
 int Verbosity(){return 0;} int FillVars();int FillCentralityInfo();float getVertexScale();
 const int NDIVS=100; const float mbd_charge_cut=.5f,mbd_time_cut=25;
 float m_mbd_total_charge=0; double m_centrality_scale=1;
 std::array<float,100> m_centrality_map;
 std::vector<std::pair<std::pair<float,float>,float>> m_vertex_scales;
 MbdPmtContainer* m_mbd_container; MbdPmtHit* m_mbd_hit; MbdOut* m_mbd_out;CentralityInfo* m_central;
};
'''
        code += "\n".join(method(sig) for sig in ("int CentralityReco::FillVars()", "int CentralityReco::FillCentralityInfo()", "float CentralityReco::getVertexScale()"))
        code += '''
uint32_t bits(float x){uint32_t u;std::memcpy(&u,&x,4);return u;}
int main(){MbdPmtContainer p;MbdOut o;CentralityInfo c;CentralityReco r;
r.m_mbd_container=&p;r.m_mbd_out=&o;r.m_central=&c;
double scale;float lo,hi,v;
while(std::cin>>o.z>>scale>>lo>>hi>>v){r.m_centrality_scale=scale;r.m_vertex_scales={{{lo,hi},v}};
for(auto& d:r.m_centrality_map)std::cin>>d;
for(auto& h:p.pmts)std::cin>>h.q>>h.t;
r.m_mbd_total_charge=0;r.FillVars();r.FillCentralityInfo();
std::cout<<bits(r.m_mbd_total_charge)<<" "<<c.b<<" "<<bits(c.c)<<"\\n";}}
'''
        rng = np.random.default_rng(42); cases = []
        for i in range(256):
            q = rng.uniform(-1, 100, 128).astype("float32")
            t = rng.uniform(-30, 30, 128).astype("float32")
            if i < 8:
                q[:] = [.49999997, .5, .50000006, 1][i % 4]; t[:] = 25 if i < 4 else -25
            c = dict(divisions=np.linspace(10000, -1, 100, dtype="float32"),
                     run_scale=float(rng.uniform(.6, 1.4)), vertex_ranges=[[-60., 60., float(np.float32(rng.uniform(.6, 1.4)))]])
            z = 60 if i % 2 else 0
            r = replay.calculate(q, t, [1] * 128, z, 1, c)
            line = [z, c["run_scale"], *c["vertex_ranges"][0], *c["divisions"]]
            line.extend(v for pair in zip(q, t) for v in pair)
            cases.append((" ".join(map(str, line)), r))
        with tempfile.TemporaryDirectory() as tmp:
            cpp = Path(tmp) / "native.cc"; cpp.write_text(code); binary = Path(tmp) / "native"
            subprocess.run(["clang++", "-std=c++17", "-O2", "-ffp-contract=off", str(cpp), "-o", str(binary)], check=True, capture_output=True)
            output = subprocess.check_output([str(binary)], input="\n".join(x[0] for x in cases) + "\n", text=True)
        self.assertEqual(len(output.splitlines()), len(cases))
        def bits(f):
            return struct.unpack("I", struct.pack("f", f))[0]
        for row, (_, r) in zip(output.splitlines(), cases):
            self.assertEqual(tuple(map(int, row.split())), (bits(r["scaled_charge"]), r["bin"], bits(r["centile"])))


if __name__ == "__main__":
    unittest.main()
