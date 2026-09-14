"""The producer sources in TreeProduction/ stay free of study numbers in comments
and of one user's absolute paths.  tools/scrub_producer_source.py --check is the
single definition of that rule; this test only runs it."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ProducerScrubTest(unittest.TestCase):
    def test_check_passes(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(ROOT / "tools" / "scrub_producer_source.py"), "--check"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

    def test_apply_is_idempotent(self) -> None:
        """A second application must report zero edits."""

        first = subprocess.run(
            [sys.executable, str(ROOT / "tools" / "scrub_producer_source.py")],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertIn("0 file edit(s)", first.stdout)


if __name__ == "__main__":
    unittest.main()
