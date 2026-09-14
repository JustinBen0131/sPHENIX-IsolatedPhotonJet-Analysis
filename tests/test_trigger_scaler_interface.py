"""Transport tests; no fixture result is a luminosity or native-physics acceptance."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
import uproot
import trigger_scaler_interface as interface


def fixture():
    values = {}
    for field, dtype in interface.SCHEMA.items():
        if field == "source_file_index":
            continue
        if dtype == "string":
            values[field] = np.asarray(["a" * 64] * 2, dtype=object)
        elif dtype == "64 * uint64":
            values[field] = np.arange(128, dtype="uint64").reshape(2, 64) + np.uint64(2**63 + 123)
        else:
            values[field] = np.zeros(2, dtype=dtype)
    values["snapshot_id"][:] = [1, 2]
    values["run"][:] = 69577
    values["segment"][:] = 3
    values["first_event_sequence"][:] = [1, 2]
    values["last_event_sequence"][:] = [1, 3]
    values["observed_events"][:] = [1, 2]
    values["packet_version"][:] = 3
    values["discontinuity_flags"][:] = [0, 1]
    events = {"event_id_hi": np.ones(3, dtype="uint64"),
              "event_id_lo": np.arange(3, dtype="uint64"),
              "run": np.full(3, 69577, dtype="int32"),
              "event_sequence": np.arange(1, 4, dtype="int64"),
              interface.REFERENCE: np.asarray([1, 2, 2], dtype="uint64")}
    return values, events


def write(path, values, events):
    with uproot.recreate(path) as root:
        tree = root.mktree(interface.TREE, {k: v for k, v in interface.SCHEMA.items() if k != "source_file_index"})
        if len(values["snapshot_id"]):
            tree.extend(values)
        tree = root.mktree("ReplayFoundationV1/RJEventV1", {k: str(v.dtype) for k, v in events.items()})
        tree.extend(events)


class ScalerInterfaceTests(unittest.TestCase):
    def test_lossless_roundtrip_and_direct_validation(self):
        values, events = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "source.root", Path(directory) / "output.root"
            write(source, values, events)
            with uproot.open(source) as original, uproot.recreate(output) as target:
                tree = target.mktree(interface.OUTPUT_TREE, interface.SCHEMA)
                self.assertEqual(interface.copy_scalers(original, tree, events, 0, require=True), 2)
            with uproot.open(output) as actual:
                interface.validate_copy([source], actual, require=True)
                rows = actual[interface.OUTPUT_TREE].arrays(library="np")
                self.assertEqual(rows["raw"].dtype, np.dtype("uint64"))
                np.testing.assert_array_equal(rows["raw"], values["raw"])
                np.testing.assert_array_equal(rows["discontinuity_flags"], [0, 1])
                self.assertEqual(int(rows["observed_events"].sum()), 3)

    def test_missing_or_wrong_reference_and_run_are_rejected(self):
        for defect in ("zero", "missing", "wrong_run", "wrong_sequence"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                values, events = fixture()
                source, output = Path(directory) / "source.root", Path(directory) / "output.root"
                write(source, values, events)
                if defect == "zero": events[interface.REFERENCE][0] = 0
                elif defect == "missing": events[interface.REFERENCE][0] = 99
                elif defect == "wrong_run": events["run"][0] = 68531
                else: events["event_sequence"][0] = 2
                with uproot.open(source) as original, uproot.recreate(output) as target:
                    tree = target.mktree(interface.OUTPUT_TREE, interface.SCHEMA)
                    with self.assertRaises(ValueError):
                        interface.copy_scalers(original, tree, events, 0)

    def test_invalid_packets_and_duplicate_snapshot_ids_are_rejected(self):
        for defect in ("packet", "duplicate"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                values, events = fixture()
                if defect == "packet": values["state"][1] = 1
                else: values["snapshot_id"][1] = 1
                source, output = Path(directory) / "source.root", Path(directory) / "output.root"
                write(source, values, events)
                with uproot.open(source) as original, uproot.recreate(output) as target:
                    tree = target.mktree(interface.OUTPUT_TREE, interface.SCHEMA)
                    with self.assertRaises(ValueError):
                        interface.copy_scalers(original, tree, events, 0, require=True)

    def test_corrupted_output_counter_is_detected_independently(self):
        values, events = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "source.root", Path(directory) / "output.root"
            write(source, values, events)
            values["raw"][1, 30] += np.uint64(1)
            with uproot.recreate(output) as target:
                tree = target.mktree(interface.OUTPUT_TREE, interface.SCHEMA)
                tree.extend({"source_file_index": np.zeros(2, dtype="int32"), **values})
            with uproot.open(output) as actual:
                with self.assertRaisesRegex(ValueError, "differs in raw"):
                    interface.validate_copy([source], actual)

    def test_strict_readback_rejects_invalid_packets_even_if_copied_exactly(self):
        values, events = fixture()
        values["state"][1] = 1
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "source.root", Path(directory) / "output.root"
            write(source, values, events)
            with uproot.open(source) as original, uproot.recreate(output) as target:
                tree = target.mktree(interface.OUTPUT_TREE, interface.SCHEMA)
                interface.copy_scalers(original, tree, events, 0, require=False)
            with uproot.open(output) as actual:
                interface.validate_copy([source], actual, require=False)
                with self.assertRaisesRegex(ValueError, "invalid packets"):
                    interface.validate_copy([source], actual, require=True)

    def test_strict_readback_rejects_empty_source_scalers(self):
        values, events = fixture()
        values = {key: value[:0] for key, value in values.items()}
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "source.root", Path(directory) / "output.root"
            write(source, values, events)
            with uproot.recreate(output) as target:
                target.mktree(interface.OUTPUT_TREE, interface.SCHEMA)
            with uproot.open(output) as actual:
                with self.assertRaisesRegex(ValueError, "empty"):
                    interface.validate_copy([source], actual, require=True)


if __name__ == "__main__":
    unittest.main()
