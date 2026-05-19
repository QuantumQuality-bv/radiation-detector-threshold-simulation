"""Analyze synthetic detector counts and generate threshold outputs and figures."""
from __future__ import annotations

import argparse
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import poisson

from simulate_counts import DETECTOR_ID, DURATION_MINUTES, EVENTS, RANDOM_SEED, load_detector_config

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = ROOT / "figures"
THRESHOLD_SIGMAS = [2, 3, 4, 5]


def estimate_background(background_df: pd.DataFrame) -> tuple[float, float, float]:
    """Estimate background mean, sample standard deviation, and Poisson standard deviation."""
    values = background_df["observed_counts"].astype(float)
    mean_counts = float(values.mean())
    std_counts = float(values.std(ddof=1))
    poisson_std_counts = float(np.sqrt(mean_counts))
    return mean_counts, std_counts, poisson_std_counts


def compute_thresholds(background_mean: float, sigma_counts: float) -> dict[int, float]:
    """Compute n-sigma threshold counts using the provided sigma estimate."""
    return {n: float(background_mean + n * sigma_counts) for n in THRESHOLD_SIGMAS}


def poisson_tail_probability(background_mean: float, threshold: float) -> float:
    """Return P(X > threshold) for X ~ Poisson(background_mean)."""
    # Crossing is defined as observed_counts > threshold.
    # For a non-integer threshold, floor(threshold) is the largest count that does not cross.
    return float(poisson.sf(math.floor(threshold), background_mean))


def first_crossing(counts_df: pd.DataFrame, threshold: float) -> int | None:
    """Return first crossing time in minutes for rows whose observed counts exceed a threshold."""
    crossings = counts_df.loc[counts_df["observed_counts"] > threshold, "time_min"]
    if crossings.empty:
        return None
    return int(crossings.iloc[0])


def build_threshold_crossings(
    background_df: pd.DataFrame,
    signal_df: pd.DataFrame,
    thresholds: dict[int, float],
) -> pd.DataFrame:
    """Return all threshold-crossing rows for both scenarios."""
    rows: list[dict[str, object]] = []
    combined = pd.concat([background_df, signal_df], ignore_index=True)
    for row in combined.itertuples(index=False):
        observed = float(row.observed_counts)
        for n, threshold in thresholds.items():
            if observed > threshold:
                rows.append(
                    {
                        "scenario": row.scenario,
                        "detector_id": row.detector_id,
                        "threshold_sigma": n,
                        "threshold_counts": threshold,
                        "time_min": int(row.time_min),
                        "observed_counts": int(row.observed_counts),
                        "excess_counts": observed - threshold,
                        "is_during_synthetic_event": bool(float(row.signal_expected_counts) > 0.0),
                    }
                )
    return pd.DataFrame(
        rows,
        columns=[
            "scenario",
            "detector_id",
            "threshold_sigma",
            "threshold_counts",
            "time_min",
            "observed_counts",
            "excess_counts",
            "is_during_synthetic_event",
        ],
    )


def build_tail_summary(
    background_df: pd.DataFrame,
    thresholds: dict[int, float],
    background_mean: float,
    poisson_std: float,
) -> pd.DataFrame:
    """Return analytical and empirical background threshold-crossing summary rows."""
    detector_id = str(background_df["detector_id"].iloc[0])
    background_counts = background_df["observed_counts"].astype(float)
    total_intervals = int(len(background_counts))
    rows: list[dict[str, object]] = []
    for n, threshold in thresholds.items():
        empirical_count = int((background_counts > threshold).sum())
        tail_probability = poisson_tail_probability(background_mean, threshold)
        rows.append(
            {
                "detector_id": detector_id,
                "threshold_sigma": n,
                "background_mean_counts": background_mean,
                "poisson_std_counts": poisson_std,
                "threshold_counts": threshold,
                "poisson_tail_probability": tail_probability,
                "background_total_intervals": total_intervals,
                "expected_background_crossings": total_intervals * tail_probability,
                "probability_at_least_one_background_crossing": 1.0
                - (1.0 - tail_probability) ** total_intervals,
                "empirical_background_crossing_count": empirical_count,
                "empirical_false_positive_fraction": empirical_count / total_intervals,
            }
        )
    return pd.DataFrame(rows)


def build_summary(
    background_df: pd.DataFrame,
    signal_df: pd.DataFrame,
    thresholds: dict[int, float],
    background_mean: float,
    background_std: float,
    poisson_std: float,
    tail_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Return false-positive and synthetic event crossing summary by threshold."""
    detector_id = str(background_df["detector_id"].iloc[0])
    background_counts = background_df["observed_counts"].astype(float)
    event_rows = signal_df.loc[signal_df["signal_expected_counts"] > 0].copy()
    event_counts = event_rows["observed_counts"].astype(float)
    tail_by_sigma = tail_summary.set_index("threshold_sigma")
    rows: list[dict[str, object]] = []
    for n, threshold in thresholds.items():
        background_crossing_count = int((background_counts > threshold).sum())
        signal_crossing_count = int((event_counts > threshold).sum())
        first_signal_crossing_min = first_crossing(signal_df, threshold)
        first_event_crossing_min = first_crossing(event_rows, threshold)
        tail_row = tail_by_sigma.loc[n]
        rows.append(
            {
                "detector_id": detector_id,
                "threshold_sigma": n,
                "background_mean_counts": background_mean,
                "background_std_counts": background_std,
                "poisson_std_counts": poisson_std,
                "threshold_counts": threshold,
                "background_crossing_count": background_crossing_count,
                "background_total_intervals": int(len(background_counts)),
                "false_positive_fraction": background_crossing_count / len(background_counts),
                "poisson_tail_probability": float(tail_row["poisson_tail_probability"]),
                "expected_background_crossings_per_120": float(tail_row["expected_background_crossings"]),
                "probability_at_least_one_background_crossing_per_120": float(
                    tail_row["probability_at_least_one_background_crossing"]
                ),
                "signal_crossing_count": signal_crossing_count,
                "signal_total_intervals": int(len(event_counts)),
                "signal_crossing_fraction": signal_crossing_count / len(event_counts),
                "first_signal_crossing_min": first_signal_crossing_min,
                "first_event_crossing_min": first_event_crossing_min,
            }
        )
    return pd.DataFrame(rows)


def build_monte_carlo_summaries(
    signal_df: pd.DataFrame,
    thresholds: dict[int, float],
    background_mean: float,
    n_trials: int,
    random_seed: int = RANDOM_SEED + 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Simulate repeated background-only and signal-plus-background runs."""
    detector_id = str(signal_df["detector_id"].iloc[0])
    intervals_per_trial = int(len(signal_df))
    signal_expected = signal_df["signal_expected_counts"].astype(float).to_numpy()
    time_min = signal_df["time_min"].astype(int).to_numpy()
    event_mask = signal_expected > 0
    event_times = time_min[event_mask]
    event_bin_count = int(event_mask.sum())
    rng = np.random.default_rng(random_seed)

    # Keep the repeated trials simple: same mean structure, new Poisson draws.
    background_counts = rng.poisson(background_mean, size=(n_trials, intervals_per_trial))
    signal_counts = rng.poisson(background_mean + signal_expected, size=(n_trials, intervals_per_trial))
    event_counts = signal_counts[:, event_mask]

    background_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    for n, threshold in thresholds.items():
        background_crossings = background_counts > threshold
        background_crossing_count = background_crossings.sum(axis=1)
        background_false_positive_fraction = background_crossing_count / intervals_per_trial
        analytical_tail = poisson_tail_probability(background_mean, threshold)
        analytical_at_least_one = 1.0 - (1.0 - analytical_tail) ** intervals_per_trial
        background_rows.append(
            {
                "detector_id": detector_id,
                "threshold_sigma": n,
                "n_trials": n_trials,
                "intervals_per_trial": intervals_per_trial,
                "threshold_counts": threshold,
                "mean_false_positive_fraction": float(background_false_positive_fraction.mean()),
                "median_false_positive_fraction": float(np.median(background_false_positive_fraction)),
                "p05_false_positive_fraction": float(np.quantile(background_false_positive_fraction, 0.05)),
                "p95_false_positive_fraction": float(np.quantile(background_false_positive_fraction, 0.95)),
                "mean_crossing_count_per_trial": float(background_crossing_count.mean()),
                "probability_at_least_one_crossing": float((background_crossing_count > 0).mean()),
                "analytical_poisson_tail_probability": analytical_tail,
                "analytical_probability_at_least_one_crossing": analytical_at_least_one,
            }
        )

        event_crossings = event_counts > threshold
        event_crossing_count = event_crossings.sum(axis=1)
        event_crossing_fraction = event_crossing_count / event_bin_count
        has_event_crossing = event_crossing_count > 0
        first_indices = np.argmax(event_crossings, axis=1)
        first_times = np.where(has_event_crossing, event_times[first_indices], np.nan)
        event_rows.append(
            {
                "detector_id": detector_id,
                "threshold_sigma": n,
                "n_trials": n_trials,
                "event_bin_count": event_bin_count,
                "threshold_counts": threshold,
                "mean_event_crossing_fraction": float(event_crossing_fraction.mean()),
                "median_event_crossing_fraction": float(np.median(event_crossing_fraction)),
                "p05_event_crossing_fraction": float(np.quantile(event_crossing_fraction, 0.05)),
                "p95_event_crossing_fraction": float(np.quantile(event_crossing_fraction, 0.95)),
                "probability_at_least_one_event_crossing": float(has_event_crossing.mean()),
                "mean_first_crossing_time_min": float(np.nanmean(first_times)),
                "median_first_crossing_time_min": float(np.nanmedian(first_times)),
            }
        )

    return pd.DataFrame(background_rows), pd.DataFrame(event_rows)


def configure_matplotlib() -> None:
    """Set report-readable white-background defaults for generated figures."""
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "black",
            "axes.labelcolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "text.color": "black",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "axes.titlecolor": "black",
            "axes.titleweight": "bold",
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "lines.linewidth": 1.8,
        }
    )


def add_event_spans(ax: plt.Axes, signal_df: pd.DataFrame) -> None:
    """Shade contiguous synthetic elevated-count intervals."""
    event_times = signal_df.loc[signal_df["signal_expected_counts"] > 0, "time_min"].astype(int).to_numpy()
    if event_times.size == 0:
        return
    start = int(event_times[0])
    previous = int(event_times[0])
    for current in event_times[1:]:
        current = int(current)
        if current != previous + 1:
            ax.axvspan(start, previous + 1, color="0.85", alpha=0.45, label="Synthetic elevated interval")
            start = current
        previous = current
    ax.axvspan(start, previous + 1, color="0.85", alpha=0.45, label="Synthetic elevated interval")


def finish_axes(ax: plt.Axes) -> None:
    """Apply common axes polish."""
    ax.grid(True, color="0.85", linewidth=0.8)
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)


def plot_background(background_df: pd.DataFrame, background_mean: float, thresholds: dict[int, float]) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(
        background_df["time_min"],
        background_df["observed_counts"],
        marker="o",
        markersize=3.5,
        color="#1f77b4",
        label="Observed counts",
    )
    ax.axhline(background_mean, color="0.25", linestyle="--", linewidth=1.8, label="Estimated background mean")
    ax.axhline(thresholds[3], color="#d62728", linestyle=":", linewidth=2.0, label="3-sigma threshold")
    ax.axhline(thresholds[5], color="#9467bd", linestyle="-.", linewidth=2.0, label="5-sigma threshold")
    ax.set_title("Background-Only Synthetic Detector Counts", color="black", fontweight="bold", pad=12)
    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("Observed counts per 60 s")
    finish_axes(ax)
    ax.legend(loc="best", frameon=True, facecolor="white", edgecolor="0.7")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "background_only_counts.png", dpi=300, facecolor="white")
    plt.close(fig)


def plot_signal(signal_df: pd.DataFrame, thresholds: dict[int, float]) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(
        signal_df["time_min"],
        signal_df["observed_counts"],
        marker="o",
        markersize=3.5,
        color="#1f77b4",
        label="Observed counts",
    )
    add_event_spans(ax, signal_df)
    ax.axhline(thresholds[3], color="#d62728", linestyle=":", linewidth=2.0, label="3-sigma threshold")
    ax.axhline(thresholds[5], color="#9467bd", linestyle="-.", linewidth=2.0, label="5-sigma threshold")
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc="best", frameon=True, facecolor="white", edgecolor="0.7")
    ax.set_title("Signal-Plus-Background Synthetic Detector Counts", color="black", fontweight="bold", pad=12)
    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("Observed counts per 60 s")
    finish_axes(ax)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "signal_plus_background_counts.png", dpi=300, facecolor="white")
    plt.close(fig)


def plot_zoom(signal_df: pd.DataFrame, thresholds: dict[int, float]) -> None:
    zoom = signal_df.loc[(signal_df["time_min"] >= 38) & (signal_df["time_min"] <= 62)].copy()
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(
        zoom["time_min"],
        zoom["observed_counts"],
        marker="o",
        markersize=4,
        color="#1f77b4",
        label="Observed counts",
    )
    add_event_spans(ax, zoom)
    ax.axhline(thresholds[3], color="#d62728", linestyle=":", linewidth=2.0, label="3-sigma threshold")
    ax.axhline(thresholds[5], color="#9467bd", linestyle="-.", linewidth=2.0, label="5-sigma threshold")
    for n, color in [(3, "#d62728"), (5, "#9467bd")]:
        crossing_time = first_crossing(zoom, thresholds[n])
        if crossing_time is not None:
            crossing_value = float(zoom.loc[zoom["time_min"] == crossing_time, "observed_counts"].iloc[0])
            ax.scatter(
                [crossing_time],
                [crossing_value],
                s=70,
                color=color,
                edgecolor="black",
                linewidth=0.7,
                zorder=5,
                label=f"First {n}-sigma crossing",
            )
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc="best", frameon=True, facecolor="white", edgecolor="0.7")
    ax.set_title("Zoom Around First Synthetic Elevated-Count Interval", color="black", fontweight="bold", pad=12)
    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("Observed counts per 60 s")
    finish_axes(ax)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "threshold_crossing_zoom.png", dpi=300, facecolor="white")
    plt.close(fig)


def plot_false_positive(summary_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    x = summary_df["threshold_sigma"]
    empirical = summary_df["false_positive_fraction"].astype(float)
    analytical = summary_df["poisson_tail_probability"].astype(float)
    at_least_one = summary_df["probability_at_least_one_background_crossing_per_120"].astype(float)
    positive_values = pd.concat([empirical[empirical > 0], analytical[analytical > 0], at_least_one[at_least_one > 0]])
    # Log plots cannot display zeros, so zero empirical values are drawn at a small floor
    # and labeled on the plot rather than treated as nonzero probabilities.
    visual_floor = max(float(positive_values.min()) / 5.0, 1e-8)
    empirical_for_plot = empirical.mask(empirical <= 0, visual_floor)

    ax.semilogy(x, empirical_for_plot, marker="o", markersize=6, color="#1f77b4", label="Empirical fraction (0 shown at floor)")
    ax.semilogy(x, analytical, marker="s", markersize=5, color="#d62728", label="Analytical Poisson P(X > Tn)")
    ax.semilogy(x, at_least_one, marker="^", markersize=6, color="#2ca02c", label="P(at least one crossing in 120 bins)")
    for n, y_value, observed in zip(x, empirical_for_plot, empirical):
        if observed == 0:
            ax.annotate(
                "0 observed",
                xy=(n, y_value),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=8,
                color="0.25",
            )
    ax.set_title("Background Threshold-Crossing Probability\nvs Threshold Multiplier", color="black", fontweight="bold", pad=12)
    ax.set_xlabel("Threshold multiplier n")
    ax.set_ylabel("Background crossing probability / fraction")
    ax.set_xticks(THRESHOLD_SIGMAS)
    ax.set_ylim(bottom=visual_floor / 2.0)
    finish_axes(ax)
    ax.legend(loc="best", frameon=True, facecolor="white", edgecolor="0.7")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "false_positive_rate_vs_threshold.png", dpi=300, facecolor="white")
    plt.close(fig)


def write_run_manifest(background_df: pd.DataFrame, signal_df: pd.DataFrame) -> None:
    """Write a compact manifest for reproducibility and audit checks."""
    detector = load_detector_config(DETECTOR_ID)
    integration_time_s = float(detector["integration_time_s"])
    manifest = {
        "random_seed": RANDOM_SEED,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "matplotlib_version": matplotlib.__version__,
        "duration_minutes": DURATION_MINUTES,
        "integration_time_s": integration_time_s,
        "number_of_time_bins": int(len(background_df)),
        "background_cpm": float(detector["background_cpm"]),
        "event_intervals": [
            {"start_min": int(event["start_min"]), "end_min": int(event["end_min"])} for event in EVENTS
        ],
        "event_added_counts": [float(event["added_counts_per_min"]) for event in EVENTS],
        "threshold_sigmas": THRESHOLD_SIGMAS,
        "script_run_timestamp": datetime.now(timezone.utc).isoformat(),
    }
    (OUTPUT_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze synthetic detector threshold crossings.")
    parser.add_argument(
        "--n-trials",
        type=int,
        default=0,
        help="Run repeated Monte Carlo summaries with the requested number of trials.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_trials < 0:
        raise ValueError("--n-trials must be nonnegative.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    background_df = pd.read_csv(OUTPUT_DIR / "background_only_counts.csv")
    signal_df = pd.read_csv(OUTPUT_DIR / "signal_plus_background_counts.csv")

    background_mean, background_std, poisson_std = estimate_background(background_df)
    thresholds = compute_thresholds(background_mean, poisson_std)

    crossings = build_threshold_crossings(background_df, signal_df, thresholds)
    tail_summary = build_tail_summary(background_df, thresholds, background_mean, poisson_std)
    summary = build_summary(
        background_df,
        signal_df,
        thresholds,
        background_mean,
        background_std,
        poisson_std,
        tail_summary,
    )

    crossings.to_csv(OUTPUT_DIR / "threshold_crossings.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "false_positive_summary.csv", index=False)
    tail_summary.to_csv(OUTPUT_DIR / "background_threshold_tail_summary.csv", index=False)
    write_run_manifest(background_df, signal_df)

    if args.n_trials > 0:
        background_mc, event_mc = build_monte_carlo_summaries(
            signal_df=signal_df,
            thresholds=thresholds,
            background_mean=background_mean,
            n_trials=args.n_trials,
        )
        background_mc.to_csv(OUTPUT_DIR / "monte_carlo_false_positive_summary.csv", index=False)
        event_mc.to_csv(OUTPUT_DIR / "monte_carlo_event_detection_summary.csv", index=False)

    plot_background(background_df, background_mean, thresholds)
    plot_signal(signal_df, thresholds)
    plot_zoom(signal_df, thresholds)
    plot_false_positive(summary)

    print("Threshold analysis complete.")
    print(f"Estimated background mean: {background_mean:.2f} counts per interval")
    print(f"Estimated background sample standard deviation: {background_std:.2f} counts")
    print(f"Poisson standard deviation used for thresholds: {poisson_std:.2f} counts")
    for n, threshold in thresholds.items():
        first = first_crossing(signal_df, threshold)
        first_label = "none" if first is None else f"{first} min"
        print(f"{n}-sigma threshold: {threshold:.2f} counts; first signal crossing: {first_label}")
    print(f"Rows written: {len(crossings)} threshold crossings, {len(summary)} summary rows")
    if args.n_trials > 0:
        print(f"Monte Carlo trials written: {args.n_trials}")


if __name__ == "__main__":
    main()
