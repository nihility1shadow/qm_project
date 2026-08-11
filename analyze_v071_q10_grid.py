#!/usr/bin/env python3
"""Rank v0.71 parameter candidates against their same-parameter QM runs."""

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


WINDOWS = ((0.0, 50.0), (50.0, 100.0), (100.0, 150.0), (150.0, 200.0))
ACTIVE_ORBITALS = np.asarray((0, 5, 6, 7, 8, 9), dtype=int)


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(signal: np.ndarray, error: np.ndarray) -> float:
    error_rms = rms(error)
    return rms(signal) / error_rms if error_rms > 0.0 else float("inf")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denominator) if denominator > 0.0 else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--target-q", type=float, default=10.0)
    parser.add_argument("--norb", type=int, default=10)
    parser.add_argument("--nel", type=int, default=5)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with args.candidates.open(newline="", encoding="utf-8-sig") as stream:
        candidates = list(csv.DictReader(stream))

    rows: list[dict[str, object]] = []
    for candidate in candidates:
        case_id = candidate["case_id"]
        qm_path = args.data_root / case_id / "qm" / "ahm-qm-s10-n5.dat"
        poisson_path = (
            args.data_root / case_id / "poisson"
            / f"ahm-sepmb-s10-n5-{candidate['ntraj']}.dat"
        )
        qm = load_dat(qm_path)
        simulation = load_dat(poisson_path)
        times = simulation[:, 0]
        qm_occ = interpolate_occupations(qm, times, args.norb)
        initial = qm_occ[0]
        qm_signal = qm_occ - initial
        simulation_signal = simulation[:, 4 : 4 + args.norb] - initial
        error = simulation_signal - qm_signal

        row: dict[str, object] = dict(candidate)
        row["qm_path"] = str(qm_path.resolve())
        row["poisson_path"] = str(poisson_path.resolve())
        row["Q_active_0_200"] = quality(
            qm_signal[:, ACTIVE_ORBITALS], error[:, ACTIVE_ORBITALS]
        )
        row["amplitude_ratio_active"] = (
            rms(simulation_signal[:, ACTIVE_ORBITALS])
            / rms(qm_signal[:, ACTIVE_ORBITALS])
        )
        row["cosine_active"] = cosine(
            qm_signal[:, ACTIVE_ORBITALS].ravel(),
            simulation_signal[:, ACTIVE_ORBITALS].ravel(),
        )
        row["max_particle_error"] = float(
            np.max(np.abs(np.sum(simulation[:, 4 : 4 + args.norb], axis=1) - args.nel))
        )

        first_three_q: list[float] = []
        first_three_strict_q: list[float] = []
        for start, stop in WINDOWS:
            mask = (times >= start - 1.0e-12) & (times <= stop + 1.0e-12)
            label = f"{int(start)}_{int(stop)}"
            q_window = quality(
                qm_signal[mask][:, ACTIVE_ORBITALS],
                error[mask][:, ACTIVE_ORBITALS],
            )
            orbital_q = [
                quality(qm_signal[mask, orbital], error[mask, orbital])
                for orbital in ACTIVE_ORBITALS
            ]
            row[f"Q_active_{label}"] = q_window
            row[f"min_orbital_Q_{label}"] = min(orbital_q)
            if stop <= 150.0:
                first_three_q.append(q_window)
                first_three_strict_q.append(min(orbital_q))

        min_q = min(first_three_q)
        multiplier = (args.target_q / min_q) ** 2 if min_q > 0.0 else float("inf")
        row["min_Q_active_0_150"] = min_q
        row["strict_min_orbital_Q_0_150"] = min(first_three_strict_q)
        row["target_Q"] = args.target_q
        row["projected_multiplier"] = multiplier
        row["projected_ntraj_Q10"] = float(candidate["ntraj"]) * multiplier
        row["projected_seconds_same_ranks_Q10"] = (
            float(candidate["poisson_seconds"]) * multiplier
        )
        row["passes_Q10"] = min_q > args.target_q
        rows.append(row)

    rows.sort(key=lambda item: float(item["min_Q_active_0_150"]), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank

    fieldnames = list(rows[0].keys())
    with (args.out_dir / "grid_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    (args.out_dir / "grid_metrics.json").write_text(
        json.dumps(rows, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    labels = [str(row["case_id"]) for row in rows]
    x = np.arange(len(rows))
    fig, axis = plt.subplots(figsize=(12, 6.2))
    for start, stop in WINDOWS[:3]:
        key = f"Q_active_{int(start)}_{int(stop)}"
        axis.plot(x, [float(row[key]) for row in rows], marker="o", label=f"{start:g}-{stop:g}")
    axis.axhline(args.target_q, color="black", linestyle="--", linewidth=1.0,
                 label=f"target Q={args.target_q:g}")
    axis.set_xticks(x, labels, rotation=30, ha="right")
    axis.set_ylabel("active-orbital aggregate Q")
    axis.set_yscale("log")
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False, ncol=4)
    fig.tight_layout()
    fig.savefig(args.out_dir / "window_q_ranking.png", dpi=220, facecolor="white")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9, 6.5))
    scatter = axis.scatter(
        [float(row["wc_eV"]) for row in rows],
        [float(row["eta"]) for row in rows],
        c=[float(row["min_Q_active_0_150"]) for row in rows],
        s=100,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.5,
    )
    for row in rows:
        axis.annotate(str(row["case_id"]), (float(row["wc_eV"]), float(row["eta"])),
                      xytext=(5, 4), textcoords="offset points", fontsize=8)
    axis.set_xlabel("wc (eV)")
    axis.set_ylabel("eta")
    axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    axis.grid(alpha=0.2)
    fig.colorbar(scatter, ax=axis, label="minimum Q through 150 a.u.")
    fig.tight_layout()
    fig.savefig(args.out_dir / "wc_eta_min_q.png", dpi=220, facecolor="white")
    plt.close(fig)

    for row in rows:
        print(
            f"#{row['rank']} {row['case_id']}: minQ={float(row['min_Q_active_0_150']):.4g} "
            f"Q50={float(row['Q_active_0_50']):.4g} "
            f"Q100={float(row['Q_active_50_100']):.4g} "
            f"Q150={float(row['Q_active_100_150']):.4g} "
            f"N(Q10)={float(row['projected_ntraj_Q10']):.4g} "
            f"runtime{row['mpi_ranks']}="
            f"{float(row['projected_seconds_same_ranks_Q10']) / 60.0:.2f} min"
        )


if __name__ == "__main__":
    main()
