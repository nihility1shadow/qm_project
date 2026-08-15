#!/usr/bin/env python3
"""Plot a long-time all-orbital QM/SepMB comparison and window Q values."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_short_qm_accuracy import interpolate_occupations, load_dat


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--simulation", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--tmax", default=500.0, type=float)
    parser.add_argument("--norb", default=10, type=int)
    args = parser.parse_args()

    simulation = load_dat(args.simulation)
    simulation = simulation[simulation[:, 0] <= args.tmax + 1.0e-12]
    times = simulation[:, 0]
    qm = load_dat(args.qm)
    qm_occ = interpolate_occupations(qm, times, args.norb)
    initial = qm_occ[0]
    qm_signal = qm_occ - initial
    simulation_signal = simulation[:, 4 : 4 + args.norb] - initial

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 5, figsize=(17, 8.6), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        combined = np.concatenate(
            (qm_signal[:, orbital], simulation_signal[:, orbital])
        )
        limit = max(float(np.max(np.abs(combined)))*1.08, 1.0e-12)
        axis.plot(
            times, simulation_signal[:, orbital], color="#d95f02",
            linewidth=0.85, alpha=0.72, label="SepMB", zorder=2,
        )
        axis.plot(
            times, qm_signal[:, orbital], color="black",
            linewidth=1.45, label="QM", zorder=3,
        )
        axis.axhline(0.0, color="#6b7280", linewidth=0.65, linestyle="--")
        axis.set_ylim(-limit, limit)
        axis.set_title(f"Orbital {orbital}")
        axis.grid(alpha=0.22)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        if orbital >= 5:
            axis.set_xlabel("time (a.u.)")
        if orbital in (0, 5):
            axis.set_ylabel("occupation change")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    order = [1, 0]
    fig.legend(
        [handles[index] for index in order], [labels[index] for index in order],
        loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2, frameon=False,
    )
    fig.suptitle(args.label, y=0.995, fontsize=13)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.925))
    fig.savefig(args.out_dir / "all_orbitals_qm_vs_sepmb.png", dpi=220,
                facecolor="white")
    plt.close(fig)

    with args.metrics.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    labels = [f"{float(row['window_start']):g}-{float(row['window_stop']):g}"
              for row in rows]
    aggregate = np.asarray([float(row["Q_active"]) for row in rows])
    individual = np.asarray([float(row["min_active_orbital_Q"]) for row in rows])
    x = np.arange(len(rows))
    fig, axis = plt.subplots(figsize=(11.5, 6.2))
    axis.plot(x, aggregate, color="#1b6ca8", marker="o", linewidth=1.8,
              label="aggregate active-orbital Q")
    axis.plot(x, individual, color="#d95f02", marker="s", linewidth=1.5,
              label="minimum individual active-orbital Q")
    axis.axhline(10.0, color="black", linewidth=1.0, linestyle="--", label="Q=10")
    axis.set_xticks(x, labels, rotation=30, ha="right")
    axis.set_yscale("log")
    axis.set_xlabel("time window (a.u.)")
    axis.set_ylabel("Q")
    axis.set_title(args.label)
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.out_dir / "window_q.png", dpi=220, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
