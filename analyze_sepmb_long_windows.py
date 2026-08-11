#!/usr/bin/env python3
"""Measure SepMB accuracy against QM over uniform long-time windows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_short_qm_accuracy import interpolate_occupations, load_dat


DEFAULT_ACTIVE = (0, 5, 6, 7, 8, 9)


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(signal: np.ndarray, error: np.ndarray) -> float:
    denominator = rms(error)
    return rms(signal) / denominator if denominator > 0.0 else float("inf")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denominator) if denominator > 0.0 else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--simulation", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--tmax", required=True, type=float)
    parser.add_argument("--window", default=50.0, type=float)
    parser.add_argument("--norb", default=10, type=int)
    parser.add_argument("--nel", default=5, type=int)
    parser.add_argument("--active", default=",".join(map(str, DEFAULT_ACTIVE)))
    parser.add_argument("--ntraj", type=int)
    parser.add_argument("--target-q", default=10.0, type=float)
    args = parser.parse_args()

    active = np.asarray([int(value) for value in args.active.split(",")], dtype=int)
    simulation = load_dat(args.simulation)
    simulation = simulation[simulation[:, 0] <= args.tmax + 1.0e-12]
    times = simulation[:, 0]
    qm = load_dat(args.qm)
    qm_occ = interpolate_occupations(qm, times, args.norb)
    initial = qm_occ[0]
    qm_signal = qm_occ - initial
    simulation_occ = simulation[:, 4 : 4 + args.norb]
    simulation_signal = simulation_occ - initial
    error = simulation_signal - qm_signal

    rows: list[dict[str, object]] = []
    start = 0.0
    while start < args.tmax - 1.0e-12:
        stop = min(start + args.window, args.tmax)
        mask = (times >= start - 1.0e-12) & (times <= stop + 1.0e-12)
        signal = qm_signal[mask][:, active]
        estimate = simulation_signal[mask][:, active]
        active_error = error[mask][:, active]
        orbital_q = [
            quality(qm_signal[mask, orbital], error[mask, orbital])
            for orbital in active
        ]
        signal_rms = rms(signal)
        q_active = quality(signal, active_error)
        multiplier = (args.target_q / q_active) ** 2 if q_active > 0.0 else float("inf")
        rows.append(
            {
                "window_start": start,
                "window_stop": stop,
                "n_points": int(np.count_nonzero(mask)),
                "Q_active": q_active,
                "min_active_orbital_Q": min(orbital_q),
                "signal_rms_active": signal_rms,
                "error_rms_active": rms(active_error),
                "amplitude_ratio_active": rms(estimate) / signal_rms,
                "cosine_active": cosine(signal.ravel(), estimate.ravel()),
                "max_particle_error": float(
                    np.max(np.abs(np.sum(simulation_occ[mask], axis=1) - args.nel))
                ),
                "target_Q": args.target_q,
                "trajectory_multiplier_for_target": multiplier,
                "projected_ntraj_for_target": (
                    args.ntraj * multiplier if args.ntraj is not None else ""
                ),
            }
        )
        start = stop

    summary = {
        "qm": str(args.qm.resolve()),
        "simulation": str(args.simulation.resolve()),
        "tmax": float(times[-1]),
        "dt": float(np.median(np.diff(times))),
        "n_points": int(len(times)),
        "active_orbitals": active.tolist(),
        "minimum_Q_active": min(float(row["Q_active"]) for row in rows),
        "windows": rows,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "long_window_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.out_dir / "long_window_metrics.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    labels = [f"{row['window_start']:g}-{row['window_stop']:g}" for row in rows]
    x = np.arange(len(rows))
    fig, axis = plt.subplots(figsize=(12, 6.4))
    axis.plot(x, [float(row["Q_active"]) for row in rows], marker="o",
              linewidth=1.8, label="aggregate active Q")
    axis.plot(x, [float(row["min_active_orbital_Q"]) for row in rows], marker="s",
              linewidth=1.3, label="minimum individual active-orbital Q")
    axis.axhline(args.target_q, color="black", linestyle="--", linewidth=1.0,
                 label=f"Q={args.target_q:g}")
    axis.set_xticks(x, labels, rotation=35, ha="right")
    axis.set_yscale("log")
    axis.set_xlabel("time window (a.u.)")
    axis.set_ylabel("Q")
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.out_dir / "long_window_q.png", dpi=220, facecolor="white")
    plt.close(fig)

    for row in rows:
        print(
            f"{row['window_start']:g}-{row['window_stop']:g}: "
            f"Q={float(row['Q_active']):.5g} "
            f"min_orbital_Q={float(row['min_active_orbital_Q']):.5g} "
            f"amp={float(row['amplitude_ratio_active']):.5g} "
            f"cos={float(row['cosine_active']):.5g}"
        )


if __name__ == "__main__":
    main()
