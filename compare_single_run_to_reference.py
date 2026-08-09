#!/usr/bin/env python3
"""Compare one large Monte Carlo run with repeated reference runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


def load_data(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="#")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path, nargs="+")
    parser.add_argument("--single", required=True, type=Path)
    parser.add_argument("--reference-ntraj", required=True, type=int)
    parser.add_argument("--single-ntraj", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--window-width", type=float, default=10.0)
    args = parser.parse_args()

    reference_data = [load_data(path) for path in args.reference]
    single_data = load_data(args.single)
    time = single_data[:, 0]
    for path, data in zip(args.reference, reference_data):
        if data.shape != single_data.shape or not np.allclose(data[:, 0], time):
            raise ValueError(f"time grid or shape differs in {path}")

    reference_occ = np.stack([data[:, 4:] for data in reference_data])
    single_occ = single_data[:, 4:]
    reference_delta = reference_occ - reference_occ[:, :1, :]
    single_delta = single_occ - single_occ[:1, :]
    reference_mean = reference_delta.mean(axis=0)
    reference_std = reference_delta.std(axis=0, ddof=1)
    difference = single_delta - reference_mean
    orbital_count = single_delta.shape[1]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    expected_ratio = math.sqrt(
        1.0 / len(reference_data) + args.reference_ntraj / args.single_ntraj
    )
    metrics_path = args.output_dir / "single_vs_reference_windows.csv"
    with metrics_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "window_start", "window_end", "orbital", "difference_rms",
            "reference_run_sd_rms", "reference_signal_rms", "cross_Q",
            "observed_sd_ratio", "expected_sd_ratio",
        ])
        start = float(time[0])
        while start < float(time[-1]) - 1e-12:
            end = min(start + args.window_width, float(time[-1]))
            mask = (time >= start) & (time <= end)
            for orbital in range(orbital_count):
                difference_rms = float(np.sqrt(np.mean(difference[mask, orbital] ** 2)))
                sd_rms = float(np.sqrt(np.mean(reference_std[mask, orbital] ** 2)))
                signal_rms = float(
                    np.sqrt(np.mean(reference_mean[mask, orbital] ** 2))
                )
                writer.writerow([
                    start, end, orbital, difference_rms, sd_rms, signal_rms,
                    signal_rms / max(difference_rms, 1e-300),
                    difference_rms / max(sd_rms, 1e-300), expected_ratio,
                ])
            start = end

    columns = min(2, orbital_count)
    rows = math.ceil(orbital_count / columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(12, max(7.8, 3.4 * rows + 1.2)),
        sharex=True, squeeze=False,
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.09, top=0.80,
                        hspace=0.40, wspace=0.30)
    colors = plt.get_cmap("tab10").colors
    for orbital, axis in enumerate(axes.flat[:orbital_count]):
        lower = reference_mean[:, orbital] - reference_std[:, orbital]
        upper = reference_mean[:, orbital] + reference_std[:, orbital]
        axis.fill_between(time, lower, upper, color=colors[orbital], alpha=0.18,
                          linewidth=0, label="30m mean +/- 1 run SD")
        axis.plot(time, reference_mean[:, orbital], color=colors[orbital],
                  linewidth=1.5, label="30m x 2 mean")
        axis.plot(time, single_delta[:, orbital], color="#111827", linewidth=1.1,
                  alpha=0.9, label="100m single run")
        axis.axhline(0, color="#64748B", linewidth=0.7, linestyle="--")
        axis.set_xlim(float(time[0]), float(time[-1]))
        axis.set_title(f"Orbital {orbital}")
        axis.grid(True, color="#CBD5E1", linewidth=0.55, alpha=0.65)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))
        axis.yaxis.set_major_formatter(formatter)
        if orbital >= (rows - 1) * columns:
            axis.set_xlabel("time (a.u.)")
        if orbital % columns == 0:
            axis.set_ylabel("occupation change")
        if orbital == 0:
            axis.legend(fontsize=8.5)
    for axis in axes.flat[orbital_count:]:
        axis.set_visible(False)

    fig.suptitle("100m single run versus two 30m reference runs", fontsize=16, y=0.97)
    fig.text(
        0.5, 0.925,
        f"wc=0.25 eV | eta=6e-5 | expected difference/30m-SD={expected_ratio:.3f}",
        ha="center", fontsize=10.5, color="#334155",
    )
    fig.savefig(args.output_dir / "all_orbitals_30m_mean_vs_100m.png",
                dpi=220, facecolor="white")
    plt.close(fig)

    particle_error = float(np.max(np.abs(single_occ.sum(axis=1) - single_data[:, 1])))
    print(f"particle_number_max_error={particle_error:.16e}")
    print(f"expected_sd_ratio={expected_ratio:.8f}")
    print(metrics_path)


if __name__ == "__main__":
    main()
