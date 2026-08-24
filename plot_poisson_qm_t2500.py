#!/usr/bin/env python3
"""Plot one stochastic Poisson run against its matching QM reference."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


def load_rows(path: Path) -> np.ndarray:
    rows = np.loadtxt(path, comments="#")
    if rows.ndim != 2 or rows.shape[1] < 14:
        raise ValueError(f"expected at least 14 columns in {path}, got {rows.shape}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--simulation", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--title", required=True)
    parser.add_argument("--subtitle", required=True)
    args = parser.parse_args()

    qm = load_rows(args.qm)
    simulation = load_rows(args.simulation)
    tmax = min(float(qm[-1, 0]), float(simulation[-1, 0]))
    qm = qm[qm[:, 0] <= tmax + 1e-12]
    simulation = simulation[simulation[:, 0] <= tmax + 1e-12]

    qm_delta = qm[:, 4:14] - qm[0, 4:14]
    simulation_delta = simulation[:, 4:14] - simulation[0, 4:14]
    colors = plt.get_cmap("tab10").colors

    fig, axes = plt.subplots(2, 5, figsize=(16, 8.8), sharex=True)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.09, top=0.84,
                        hspace=0.42, wspace=0.31)
    for orbital, axis in enumerate(axes.flat):
        axis.plot(qm[:, 0], qm_delta[:, orbital], color="#111827",
                  linewidth=1.35, label="QM")
        axis.plot(simulation[:, 0], simulation_delta[:, orbital],
                  color=colors[orbital], linewidth=1.15, alpha=0.9,
                  label="single Poisson run")
        axis.axhline(0.0, color="#94a3b8", linewidth=0.65, linestyle="--")
        values = np.concatenate((qm_delta[:, orbital], simulation_delta[:, orbital]))
        low = float(np.min(values))
        high = float(np.max(values))
        span = max(high - low, 1e-14)
        axis.set_xlim(0.0, tmax)
        axis.set_ylim(low - 0.08 * span, high + 0.08 * span)
        state = "occupied" if orbital < 5 else "empty"
        axis.set_title(f"Orbital {orbital} ({state})", fontsize=10)
        axis.grid(True, color="#d1d5db", linewidth=0.5, alpha=0.65)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))
        axis.yaxis.set_major_formatter(formatter)
        if orbital >= 5:
            axis.set_xlabel("time (a.u.)")
        if orbital in (0, 5):
            axis.set_ylabel("occupation change")
    axes.flat[0].legend(loc="best", fontsize=8, frameon=False)
    fig.suptitle(args.title, fontsize=16, y=0.975)
    fig.text(0.5, 0.93, args.subtitle, ha="center", fontsize=10.5,
             color="#334155")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    orbital_path = args.output_dir / "all_orbitals_poisson_vs_qm.png"
    fig.savefig(orbital_path, dpi=220, facecolor="white")
    plt.close(fig)

    with args.metrics.open(newline="", encoding="utf-8") as handle:
        metrics = list(csv.DictReader(handle))
    midpoint = np.array([(float(row["start"]) + float(row["stop"])) / 2
                         for row in metrics])
    aggregate = np.array([float(row["aggregate_Q"]) for row in metrics])
    strict = np.array([float(row["min_orbital_Q"]) for row in metrics])

    fig, axis = plt.subplots(figsize=(11.5, 5.7))
    axis.semilogy(midpoint, aggregate, marker="o", markersize=3.8,
                  linewidth=1.4, color="#2563eb", label="aggregate Q")
    axis.semilogy(midpoint, strict, marker="s", markersize=3.4,
                  linewidth=1.35, color="#dc2626", label="strict orbital Q")
    axis.axhline(10.0, color="#111827", linestyle="--", linewidth=1.0,
                 label="Q = 10")
    axis.axhline(1.0, color="#64748b", linestyle=":", linewidth=1.0,
                 label="Q = 1")
    axis.set_xlim(0.0, tmax)
    axis.set_xlabel("time-window midpoint (a.u.)")
    axis.set_ylabel("Q = RMS(QM signal) / RMS(error)")
    axis.set_title("Single-run accuracy by 100 a.u. window")
    axis.grid(True, which="both", color="#d1d5db", linewidth=0.55, alpha=0.7)
    axis.legend(frameon=False, ncol=4, loc="upper right")
    q_path = args.output_dir / "q_by_window.png"
    fig.tight_layout()
    fig.savefig(q_path, dpi=220, facecolor="white")
    plt.close(fig)

    print(orbital_path)
    print(q_path)


if __name__ == "__main__":
    main()
