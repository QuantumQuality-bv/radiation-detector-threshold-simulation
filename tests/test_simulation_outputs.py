"""Basic checks for generated synthetic detector count outputs."""
from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


class SimulationOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, "src/simulate_counts.py"], cwd=ROOT, check=True)

    def test_count_csvs_have_expected_row_counts(self) -> None:
        background = pd.read_csv(ROOT / "outputs" / "background_only_counts.csv")
        signal = pd.read_csv(ROOT / "outputs" / "signal_plus_background_counts.csv")
        self.assertEqual(len(background), 120)
        self.assertEqual(len(signal), 120)

    def test_synthetic_event_intervals_have_signal_counts(self) -> None:
        signal = pd.read_csv(ROOT / "outputs" / "signal_plus_background_counts.csv")
        event1 = signal.loc[(signal["time_min"] >= 45) & (signal["time_min"] < 55)]
        event2 = signal.loc[(signal["time_min"] >= 85) & (signal["time_min"] < 90)]
        self.assertTrue((event1["signal_expected_counts"] > 0).all())
        self.assertTrue((event2["signal_expected_counts"] > 0).all())
        self.assertEqual(len(event1), 10)
        self.assertEqual(len(event2), 5)


if __name__ == "__main__":
    unittest.main()
