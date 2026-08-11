#!/usr/bin/env python3
"""Compare QM, continuous SepMB, and sequential-window SepMB results."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_dat(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    return data[None, :] if data.ndim == 1 else data


def interpolate_occupations(data: np.ndarray, times: np.ndarray, norb: int) -> np.ndarray:
    result = np.empty((len(times), norb))
    for orbital in range(norb):
        result[:, orbital] = np.interp(times, data[:, 0], data[:, 4 + orbital])
    return result


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(signal: np.ndarray, error: np.ndarray) -> float:
    error_rms = rms(error)
    return rms(signal) / error_rms if error_rms > 0.0 else float("inf")


def parse_timing(path: Path) -> tuple[float, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    real_match = re.search(r"^real\s+([0-9.]+)$", text, re.MULTILINE)
    user_match = re.search(r"^user\s+([0-9.]+)$", text, re.MULTILINE)
    return (
        float(real_match.group(1)) if real_match else float("nan"),
        float(user_match.group(1)) if user_match else float("nan"),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--continuous", required=True, type=Path)
    parser.add_argument("--windowed", required=True, type=Path)
    parser.add_argument("--qm-log", required=True, type=Path)
    parser.add_argument("--continuous-log", required=True, type=Path)
    parser.add_argument("--windowed-log", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--norb", default=10, type=int)
    parser.add_argument("--nel", default=5, type=int)
    parser.add_argument("--ntraj", default=100000, type=int)
    parser.add_argument("--window", default=25.0, type=float)
    parser.add_argument("--tmax", default=100.0, type=float)
    args = parser.parse_args()

    qm_data = load_dat(args.qm)
    continuous_data = load_dat(args.continuous)
    windowed_data = load_dat(args.windowed)
    continuous_data = continuous_data[continuous_data[:, 0] <= args.tmax + 1.0e-12]
    windowed_data = windowed_data[windowed_data[:, 0] <= args.tmax + 1.0e-12]
    times = continuous_data[:, 0]
    if not np.allclose(times, windowed_data[:, 0]):
        raise ValueError("continuous and windowed time grids differ")

    qm_occ = interpolate_occupations(qm_data, times, args.norb)
    continuous_occ = continuous_data[:, 4 : 4 + args.norb]
    windowed_occ = windowed_data[:, 4 : 4 + args.norb]
    initial = qm_occ[0]
    qm_signal = qm_occ - initial
    continuous_signal = continuous_occ - initial
    windowed_signal = windowed_occ - initial

    rows: list[dict[str, float | int | str]] = []
    start = 0.0
    while start < args.tmax - 1.0e-12:
        stop = min(start + args.window, args.tmax)
        mask = (times >= start - 1.0e-12) & (times <= stop + 1.0e-12)
        reference = qm_signal[mask]
        for method, estimate in (
            ("continuous_v071", continuous_signal[mask]),
            ("windowed_v075", windowed_signal[mask]),
        ):
            q_value = quality(reference, estimate - reference)
            orbital_q = [
                quality(reference[:, orbital], estimate[:, orbital] - reference[:, orbital])
                for orbital in range(args.norb)
            ]
            rows.append(
                {
                    "method": method,
                    "window_start": start,
                    "window_stop": stop,
                    "Q_all_orbitals": q_value,
                    "minimum_individual_orbital_Q": min(orbital_q),
                    "signal_rms": rms(reference),
                    "error_rms": rms(estimate - reference),
                    "amplitude_ratio": rms(estimate) / rms(reference),
                    "projected_ntraj_for_Q10": args.ntraj * (10.0 / q_value) ** 2,
                    "max_particle_error": float(
                        np.max(np.abs(np.sum(estimate + initial, axis=1) - args.nel))
                    ),
                }
            )
        start = stop

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "window_comparison_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    qm_real, qm_user = parse_timing(args.qm_log)
    continuous_real, continuous_user = parse_timing(args.continuous_log)
    windowed_real, windowed_user = parse_timing(args.windowed_log)
    timing = {
        "qm": {"wall_seconds": qm_real, "user_cpu_seconds": qm_user, "mpi_ranks": 1},
        "continuous_v071": {
            "wall_seconds": continuous_real,
            "user_cpu_seconds": continuous_user,
            "mpi_ranks": 128,
            "forward_trajectories": args.ntraj,
            "back_replicas": 64,
        },
        "windowed_v075": {
            "wall_seconds": windowed_real,
            "user_cpu_seconds": windowed_user,
            "mpi_ranks": 128,
            "ket_trajectories": args.ntraj,
            "bra_trajectories": args.ntraj,
            "pair_replicas": 8,
            "window_steps": 50,
            "window_au": 25.0,
        },
    }
    summary = {
        "parameters": {
            "wc_eV": 3.389460,
            "eta": 2.961e-5,
            "delE_eV": -17.194874968839816,
            "dt_au": 0.5,
            "tmax_au": args.tmax,
            "norb": args.norb,
            "nel": args.nel,
        },
        "timing": timing,
        "metrics": rows,
    }
    (args.out_dir / "window_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    colors = {"QM": "#111111", "continuous": "#1874b4", "windowed": "#d1495b"}
    fig, axes = plt.subplots(5, 2, figsize=(13.5, 15.5), sharex=True)
    for orbital, axis in enumerate(axes.ravel()):
        axis.plot(times, qm_signal[:, orbital], color=colors["QM"], linewidth=1.8, label="QM")
        axis.plot(times, continuous_signal[:, orbital], color=colors["continuous"],
                  linewidth=1.25, label="continuous v0.71")
        axis.plot(times, windowed_signal[:, orbital], color=colors["windowed"],
                  linewidth=1.05, label="4 x 25 a.u. v0.75")
        for boundary in (25.0, 50.0, 75.0):
            axis.axvline(boundary, color="#777777", linestyle=":", linewidth=0.8)
        axis.set_title(f"orbital {orbital}")
        axis.set_ylabel(r"$\Delta n_j(t)$")
        axis.grid(alpha=0.2)
    for axis in axes[-1]:
        axis.set_xlabel("time (a.u.)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("100 a.u. comparison, 100,000 trajectories", y=0.997)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.981),
               ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.952))
    fig.savefig(args.out_dir / "all_orbitals_qm_continuous_windowed.png", dpi=220,
                facecolor="white")
    plt.close(fig)

    labels = ["0-25", "25-50", "50-75", "75-100"]
    continuous_q = [float(row["Q_all_orbitals"]) for row in rows
                    if row["method"] == "continuous_v071"]
    windowed_q = [float(row["Q_all_orbitals"]) for row in rows
                  if row["method"] == "windowed_v075"]
    x = np.arange(len(labels))
    fig, axis = plt.subplots(figsize=(10.5, 5.8))
    axis.plot(x, continuous_q, "o-", color=colors["continuous"], linewidth=1.8,
              label="continuous v0.71")
    axis.plot(x, windowed_q, "s-", color=colors["windowed"], linewidth=1.8,
              label="4 x 25 a.u. v0.75")
    axis.axhline(10.0, color="#111111", linestyle="--", linewidth=1.0, label="Q=10")
    axis.set_xticks(x, labels)
    axis.set_yscale("log")
    axis.set_xlabel("time window (a.u.)")
    axis.set_ylabel("Q over all 10 orbitals")
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.out_dir / "window_q_comparison.png", dpi=220, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
