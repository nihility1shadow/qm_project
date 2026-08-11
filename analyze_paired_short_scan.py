#!/usr/bin/env python3
"""Rank paired SepMB/QM short-time parameter cases."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_short_qm_accuracy import (
    correlation,
    cosine_similarity,
    interpolate_occupations,
    load_dat,
    quality,
    rms,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--norb", type=int, default=10)
    parser.add_argument("--nel", type=int, default=5)
    parser.add_argument("--tmax", type=float, default=25.0)
    parser.add_argument("--rank-by", choices=("all", "min-half"), default="all")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with args.candidates.open(newline="", encoding="utf-8") as stream:
        candidates = list(csv.DictReader(stream))

    summary_rows: list[dict[str, object]] = []
    orbital_rows: list[dict[str, object]] = []
    case_data: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for candidate in candidates:
        case_id = candidate["case_id"]
        ntraj = int(candidate["ntraj"])
        poisson_path = (
            args.run_root / case_id / "poisson" /
            f"ahm-sepmb-s{args.norb}-n{args.nel}-{ntraj}.dat"
        )
        qm_path = args.run_root / case_id / "qm" / f"ahm-qm-s{args.norb}-n{args.nel}.dat"
        poisson = load_dat(poisson_path)
        qm = load_dat(qm_path)
        poisson = poisson[poisson[:, 0] <= args.tmax + 1.0e-12]
        qm = qm[qm[:, 0] <= args.tmax + 1.0e-12]
        times = poisson[:, 0]
        qm_occ = interpolate_occupations(qm, times, args.norb)
        sim_occ = poisson[:, 4 : 4 + args.norb]
        initial = qm_occ[0]
        qm_signal = qm_occ - initial
        sim_signal = sim_occ - initial
        error = sim_signal - qm_signal

        active = np.array([0, *range(args.nel, args.norb)], dtype=int)
        static = np.arange(1, args.nel, dtype=int)
        signal_rms = rms(qm_signal)
        error_rms = rms(error)
        sim_rms = rms(sim_signal)
        active_q = quality(qm_signal[:, active], error[:, active])
        static_leakage = rms(sim_signal[:, static]) if len(static) else 0.0
        half_time = 0.5 * args.tmax
        first_half = times <= half_time + 1.0e-12
        second_half = times >= half_time - 1.0e-12
        q_first_half = quality(qm_signal[first_half], error[first_half])
        q_second_half = quality(qm_signal[second_half], error[second_half])
        row: dict[str, object] = {
            **candidate,
            "Q_qm_all": quality(qm_signal, error),
            "Q_qm_active": active_q,
            "Q_qm_first_half": q_first_half,
            "Q_qm_second_half": q_second_half,
            "Q_qm_min_half": min(q_first_half, q_second_half),
            "amplitude_ratio": sim_rms / signal_rms if signal_rms > 0.0 else float("nan"),
            "cosine_similarity": cosine_similarity(qm_signal.ravel(), sim_signal.ravel()),
            "qm_signal_rms": signal_rms,
            "error_rms": error_rms,
            "static_orbital_leakage_rms": static_leakage,
            "max_abs_occupation_error": float(np.max(np.abs(error))),
            "max_particle_error": float(np.max(np.abs(np.sum(sim_occ, axis=1) - args.nel))),
        }
        summary_rows.append(row)
        case_data[case_id] = (times, qm_signal, sim_signal)

        for orbital in range(args.norb):
            q = qm_signal[:, orbital]
            s = sim_signal[:, orbital]
            e = s - q
            q_rms = rms(q)
            s_rms = rms(s)
            orbital_rows.append(
                {
                    "case_id": case_id,
                    "orbital": orbital,
                    "qm_signal_rms": q_rms,
                    "simulation_signal_rms": s_rms,
                    "error_rms": rms(e),
                    "Q_qm": quality(q, e),
                    "amplitude_ratio": s_rms / q_rms if q_rms > 0.0 else float("nan"),
                    "cosine_similarity": cosine_similarity(q, s),
                    "pearson_correlation": correlation(q, s),
                    "max_abs_error": float(np.max(np.abs(e))),
                }
            )

    rank_key = "Q_qm_all" if args.rank_by == "all" else "Q_qm_min_half"
    summary_rows.sort(key=lambda item: float(item[rank_key]), reverse=True)
    for rank, row in enumerate(summary_rows, start=1):
        row["rank"] = rank

    with (args.out_dir / "paired_qm_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    with (args.out_dir / "paired_qm_per_orbital.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(orbital_rows[0].keys()))
        writer.writeheader()
        writer.writerows(orbital_rows)

    payload = {"tmax": args.tmax, "ranked_cases": summary_rows, "orbitals": orbital_rows}
    (args.out_dir / "paired_qm_metrics.json").write_text(
        json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    labels = [str(row["case_id"]) for row in summary_rows]
    q_values = [float(row[rank_key]) for row in summary_rows]
    fig, axis = plt.subplots(figsize=(10.5, 5.8))
    bars = axis.bar(labels, q_values, color="#2f6f8f")
    axis.axhline(10.0, color="#b33a3a", linestyle="--", linewidth=1.2, label="Q = 10")
    axis.set_ylabel("Q against same-parameter QM")
    axis.set_xlabel("parameter case")
    rank_text = "whole interval" if args.rank_by == "all" else "weaker half interval"
    axis.set_title(
        f"v0.64 parameter ranking by {rank_text}, t <= {args.tmax:g} a.u."
    )
    axis.grid(axis="y", alpha=0.25)
    axis.legend(frameon=False)
    for bar, value in zip(bars, q_values):
        axis.text(bar.get_x() + bar.get_width() / 2.0, value, f"{value:.2f}",
                  ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "paired_qm_ranking.png", dpi=220, facecolor="white")
    plt.close(fig)

    best = summary_rows[0]
    best_id = str(best["case_id"])
    times, qm_signal, sim_signal = case_data[best_id]
    fig, axes = plt.subplots(2, 5, figsize=(17, 8.5), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        axis.plot(times, qm_signal[:, orbital], color="black", linewidth=2.0, label="QM")
        axis.plot(times, sim_signal[:, orbital], color="#c44e52", linewidth=1.5,
                  label="SepMB v0.64")
        axis.axhline(0.0, color="#6b7280", linewidth=0.7, linestyle="--")
        axis.set_title(f"Orbital {orbital}")
        axis.grid(alpha=0.25)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        if orbital >= 5:
            axis.set_xlabel("time (a.u.)")
        if orbital in (0, 5):
            axis.set_ylabel("occupation change")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.955))
    fig.suptitle(
        f"Best paired case {best_id}: wc={float(best['wc_eV']):.6g} eV, "
        f"eta={float(best['eta']):.6g}, Ntraj={int(best['ntraj']):,}",
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.91))
    fig.savefig(args.out_dir / "best_case_all_orbitals_vs_qm.png", dpi=220,
                facecolor="white")
    plt.close(fig)

    for row in summary_rows:
        print(
            f"{int(row['rank']):2d} {row['case_id']}: "
            f"wc={float(row['wc_eV']):.6g} eta={float(row['eta']):.6g} "
            f"Q={float(row['Q_qm_all']):.4f} "
            f"Qhalf={float(row['Q_qm_min_half']):.4f} "
            f"amp={float(row['amplitude_ratio']):.4f} "
            f"cos={float(row['cosine_similarity']):.4f}"
        )


if __name__ == "__main__":
    main()
