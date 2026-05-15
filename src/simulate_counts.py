"""Generate synthetic detector count time series for the threshold demo."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "outputs"

DETECTOR_ID = "D001"
DURATION_MINUTES = 120
RANDOM_SEED = 42

# Half-open intervals: start_min <= t < end_min. This matches how Python slices work.
EVENTS = [
    {"start_min": 45, "end_min": 55, "added_counts_per_min": 350.0},
    {"start_min": 85, "end_min": 90, "added_counts_per_min": 650.0},
]


def load_detector_config(detector_id: str = DETECTOR_ID) -> pd.Series:
    """Load one detector row from the synthetic configuration file."""
    config_path = DATA_DIR / "synthetic_detector_config.csv"
    config = pd.read_csv(config_path)
    matches = config.loc[config["detector_id"] == detector_id]
    if matches.empty:
        raise ValueError(f"Detector {detector_id!r} was not found in {config_path}.")
    return matches.iloc[0]


def signal_profile(time_min: np.ndarray, integration_time_s: float) -> np.ndarray:
    """Return expected added signal counts for each integration interval."""
    signal = np.zeros_like(time_min, dtype=float)
    integration_scale = integration_time_s / 60.0

    # Signal values are means for Poisson sampling, not deterministic observed counts.
    for event in EVENTS:
        in_event = (time_min >= event["start_min"]) & (time_min < event["end_min"])
        signal[in_event] += event["added_counts_per_min"] * integration_scale
    return signal


def build_scenario_frame(
    time_min: np.ndarray,
    detector_id: str,
    background_expected: float,
    signal_expected: np.ndarray,
    observed_counts: np.ndarray,
    scenario: str,
) -> pd.DataFrame:
    """Assemble one scenario output table."""
    return pd.DataFrame(
        {
            "time_min": time_min.astype(int),
            "detector_id": detector_id,
            "background_expected_counts": background_expected,
            "signal_expected_counts": signal_expected,
            "total_expected_counts": background_expected + signal_expected,
            "observed_counts": observed_counts.astype(int),
            "scenario": scenario,
        }
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    detector = load_detector_config()
    detector_id = str(detector["detector_id"])
    integration_time_s = float(detector["integration_time_s"])
    background_cpm = float(detector["background_cpm"])
    background_expected = background_cpm * integration_time_s / 60.0

    time_min = np.arange(DURATION_MINUTES, dtype=int)
    signal_expected = signal_profile(time_min, integration_time_s)

    # One generator is used for both sequences so the run is exactly reproducible.
    rng = np.random.default_rng(RANDOM_SEED)
    background_only_observed = rng.poisson(background_expected, size=DURATION_MINUTES)
    signal_plus_background_observed = rng.poisson(
        background_expected + signal_expected,
        size=DURATION_MINUTES,
    )

    background_only = build_scenario_frame(
        time_min=time_min,
        detector_id=detector_id,
        background_expected=background_expected,
        signal_expected=np.zeros_like(time_min, dtype=float),
        observed_counts=background_only_observed,
        scenario="background_only",
    )
    signal_plus_background = build_scenario_frame(
        time_min=time_min,
        detector_id=detector_id,
        background_expected=background_expected,
        signal_expected=signal_expected,
        observed_counts=signal_plus_background_observed,
        scenario="signal_plus_background",
    )

    background_path = OUTPUT_DIR / "background_only_counts.csv"
    signal_path = OUTPUT_DIR / "signal_plus_background_counts.csv"
    background_only.to_csv(background_path, index=False)
    signal_plus_background.to_csv(signal_path, index=False)

    print("Synthetic detector count simulation complete.")
    print(f"Rows written: {len(background_only)} background-only, {len(signal_plus_background)} signal-plus-background")
    print(f"Detector: {detector_id}")
    print(f"Background mean used: {background_expected:.2f} counts per {integration_time_s:.0f} s interval")
    for index, event in enumerate(EVENTS, start=1):
        print(
            "Event "
            f"{index}: minutes {event['start_min']} to {event['end_min']}, "
            f"added signal {event['added_counts_per_min']:.0f} counts/min"
        )


if __name__ == "__main__":
    main()
