#!/usr/bin/env python3
"""Compare strict-Q windows at two trajectory counts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_first(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return next(csv.DictReader(handle))


def q_windows(row: dict[str, str]) -> list[tuple[float, float, float]]:
    windows = []
    for key, value in row.items():
        if not key.startswith("min_active_Q_") or not value:
            continue
        start, end = map(float, key.removeprefix("min_active_Q_").split("_"))
        windows.append((start, end, float(value)))
    return sorted(windows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-ranking", required=True, type=Path)
    parser.add_argument("--large-ranking", required=True, type=Path)
    parser.add_argument("--baseline-ntraj", required=True, type=int)
    parser.add_argument("--large-ntraj", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threshold", type=float, default=10.0)
    args = parser.parse_args()

    baseline = q_windows(load_first(args.baseline_ranking))
    large = q_windows(load_first(args.large_ranking))
    if [(a, b) for a, b, _ in baseline] != [(a, b) for a, b, _ in large]:
        raise ValueError("Q-window grids differ")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scale = np.sqrt(args.large_ntraj / args.baseline_ntraj)
    midpoint = np.array([(start + end) / 2 for start, end, _ in baseline])
    baseline_q = np.array([q for _, _, q in baseline])
    large_q = np.array([q for _, _, q in large])

    with (args.output_dir / "q_scaling_by_window.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "window_start", "window_end", "baseline_Q", "large_Q",
            "observed_ratio", "sqrt_ntraj_ratio",
        ])
        for (start, end, _), q0, q1 in zip(baseline, baseline_q, large_q):
            writer.writerow([start, end, q0, q1, q1 / q0, scale])

    fig, axis = plt.subplots(figsize=(11.5, 6.5))
    axis.plot(midpoint, baseline_q, "o-", linewidth=1.8,
              label=f"{args.baseline_ntraj:,} trajectories (measured)")
    axis.plot(midpoint, baseline_q * scale, "o--", linewidth=1.6,
              label=f"sqrt(N) prediction x {scale:.3f}")
    axis.plot(midpoint, large_q, "o-", linewidth=2.0,
              label=f"{args.large_ntraj:,} trajectories (measured)")
    axis.axhline(args.threshold, color="#B91C1C", linestyle="--", linewidth=1.4,
                 label=f"target Q > {args.threshold:g}")
    axis.set(
        xlabel="time-window midpoint (a.u.)",
        ylabel="minimum active-orbital Q",
        title="Trajectory-count scaling of strict convergence Q",
    )
    axis.set_xlim(0, max(end for _, end, _ in large))
    axis.set_yscale("log")
    axis.grid(True, color="#CBD5E1", linewidth=0.6, alpha=0.75)
    axis.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output_dir / "q_1m_vs_30m.png", dpi=220, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
