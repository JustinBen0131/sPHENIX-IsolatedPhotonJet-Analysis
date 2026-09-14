#!/usr/bin/env python3
"""Run the test suite: python tests/run_tests.py [pattern]

Stage directories are put on the import path so tests can import the stage
modules directly (photon_selection, make_histograms, photonjet, ...).
"""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
for directory in ("FinalAnalysis", "PhotonID", "TreeToHists", "TreeProduction", "TreeProduction/producer", "tests"):
    sys.path.insert(0, str(ROOT / directory))

pattern = sys.argv[1] if len(sys.argv) > 1 else "test_*.py"
suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern=pattern, top_level_dir=str(ROOT / "tests"))
result = unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(0 if result.wasSuccessful() else 1)
