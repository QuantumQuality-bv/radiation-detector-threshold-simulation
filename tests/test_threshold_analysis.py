"""Checks for threshold calculations, probability summaries, and figure artifacts."""
from __future__ import annotations

import math
import subprocess
import sys
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


class ThresholdAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, "src/simulate_counts.py"], cwd=ROOT, check=True)
        subprocess.run([sys.executable, "src/threshold_analysis.py"], cwd=ROOT, check=True)

    def test_threshold_formula_uses_poisson_standard_deviation(self) -> None:
        summary = pd.read_csv(ROOT / "outputs" / "false_positive_summary.csv")
        for row in summary.itertuples(index=False):
            expected_threshold = row.background_mean_counts + row.threshold_sigma * math.sqrt(row.background_mean_counts)
            self.assertAlmostEqual(row.threshold_counts, expected_threshold, places=9)

    def test_probability_columns_are_bounded(self) -> None:
        summary = pd.read_csv(ROOT / "outputs" / "false_positive_summary.csv")
        tail = pd.read_csv(ROOT / "outputs" / "background_threshold_tail_summary.csv")
        self.assertTrue(summary["false_positive_fraction"].between(0, 1).all())
        self.assertTrue(summary["poisson_tail_probability"].between(0, 1).all())
        self.assertTrue(summary["probability_at_least_one_background_crossing_per_120"].between(0, 1).all())
        self.assertTrue(tail["poisson_tail_probability"].between(0, 1).all())
        self.assertTrue(tail["probability_at_least_one_background_crossing"].between(0, 1).all())

    def test_default_seed_has_3sigma_synthetic_event_crossing(self) -> None:
        crossings = pd.read_csv(ROOT / "outputs" / "threshold_crossings.csv")
        event_crossings = crossings.loc[
            (crossings["threshold_sigma"] == 3) & (crossings["is_during_synthetic_event"])
        ]
        self.assertGreaterEqual(len(event_crossings), 1)

    def test_required_figure_files_exist(self) -> None:
        figure_names = [
            "background_only_counts.png",
            "signal_plus_background_counts.png",
            "threshold_crossing_zoom.png",
            "false_positive_rate_vs_threshold.png",
        ]
        image_names = [
            "background_only_counts_matlab.png",
            "signal_plus_background_counts_matlab.png",
            "threshold_crossing_zoom_matlab.png",
            "false_positive_rate_vs_threshold_matlab.png",
        ]
        for name in figure_names:
            path = ROOT / "figures" / name
            self.assertTrue(path.exists(), f"Missing figure file: {path}")
            self.assertGreater(path.stat().st_size, 0)
        for name in image_names:
            path = ROOT / "images" / name
            self.assertTrue(path.exists(), f"Missing Overleaf image file: {path}")
            self.assertGreater(path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
