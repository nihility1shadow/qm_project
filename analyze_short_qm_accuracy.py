#!/usr/bin/env python3
"""Compare short-time SepMB runs with a same-parameter QM reference."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_dat(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    return data


def interpolate_occupations(data: np.ndarray, times: np.ndarray, norb: int) -> np.ndarray:
    result = np.empty((len(times), norb), dtype=float)
    for orbital in range(norb):
        result[:, orbital] = np.interp(times, data[:, 0], data[:, 4 + orbital])
    return result


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(reference: np.ndarray, error: np.ndarray) -> float:
    denominator = rms(error)
    if denominator == 0.0:
        return float("inf")
    return rms(reference) / denominator


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denominator) if denominator > 0.0 else float("nan")


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def parse_run(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("run must use LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label:
        raise argparse.ArgumentTypeError("run label cannot be empty")
    return label, Path(raw_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--run", required=True, action="append", type=parse_run)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--norb", type=int, default=10)
    parser.add_argument("--nel", type=int, default=5)
    parser.add_argument("--tmax", type=float, default=25.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    qm = load_dat(args.qm)
    qm = qm[qm[:, 0] <= args.tmax + 1.0e-12]
    initial = qm[0, 4 : 4 + args.norb]

    summaries: list[dict[str, object]] = []
    orbital_rows: list[dict[str, object]] = []
    plot_data: list[tuple[str, np.ndarray, np.ndarray]] = []

    for label, path in args.run:
        simulation = load_dat(path)
        simulation = simulation[simulation[:, 0] <= args.tmax + 1.0e-12]
        times = simulation[:, 0]
        qm_occ = interpolate_occupations(qm, times, args.norb)
        sim_occ = simulation[:, 4 : 4 + args.norb]
        qm_signal = qm_occ - initial
        sim_signal = sim_occ - initial
        error = sim_signal - qm_signal

        signal_rms = rms(qm_signal)
        error_rms = rms(error)
        sim_rms = rms(sim_signal)
        half_time = 0.5 * args.tmax
        first_half = times <= half_time + 1.0e-12
        second_half = times >= half_time - 1.0e-12
        q_first_half = quality(qm_signal[first_half], error[first_half])
        q_second_half = quality(qm_signal[second_half], error[second_half])
        summary: dict[str, object] = {
            "label": label,
            "path": str(path.resolve()),
            "n_points": int(len(times)),
            "tmax": float(times[-1]),
            "qm_signal_rms": signal_rms,
            "simulation_signal_rms": sim_rms,
            "error_rms": error_rms,
            "Q_qm": quality(qm_signal, error),
            "Q_qm_first_half": q_first_half,
            "Q_qm_second_half": q_second_half,
            "Q_qm_min_half": min(q_first_half, q_second_half),
            "amplitude_ratio": sim_rms / signal_rms if signal_rms > 0.0 else float("nan"),
            "cosine_similarity": cosine_similarity(qm_signal.ravel(), sim_signal.ravel()),
            "max_abs_occupation_error": float(np.max(np.abs(error))),
            "max_particle_error": float(np.max(np.abs(np.sum(sim_occ, axis=1) - args.nel))),
        }
        summaries.append(summary)
        plot_data.append((label, times, sim_signal))

        for orbital in range(args.norb):
            q = qm_signal[:, orbital]
            s = sim_signal[:, orbital]
            e = s - q
            q_rms = rms(q)
            s_rms = rms(s)
            orbital_rows.append(
                {
                    "label": label,
                    "orbital": orbital,
                    "initially_occupied": orbital < args.nel,
                    "qm_signal_rms": q_rms,
                    "simulation_signal_rms": s_rms,
                    "error_rms": rms(e),
                    "Q_qm": quality(q, e),
                    "amplitude_ratio": s_rms / q_rms if q_rms > 0.0 else float("nan"),
                    "cosine_similarity": cosine_similarity(q, s),
                    "pearson_correlation": correlation(q, s),
                    "qm_peak_abs": float(np.max(np.abs(q))),
                    "simulation_peak_abs": float(np.max(np.abs(s))),
                    "max_abs_error": float(np.max(np.abs(e))),
                }
            )

    payload = {
        "qm": str(args.qm.resolve()),
        "norb": args.norb,
        "nel": args.nel,
        "tmax": args.tmax,
        "runs": summaries,
        "orbitals": orbital_rows,
    }
    (args.out_dir / "short_qm_accuracy.json").write_text(
        json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    with (args.out_dir / "short_qm_accuracy.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)

    with (args.out_dir / "short_qm_accuracy_per_orbital.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(orbital_rows[0].keys()))
        writer.writeheader()
        writer.writerows(orbital_rows)

    qm_times = plot_data[0][1]
    qm_plot_occ = interpolate_occupations(qm, qm_times, args.norb)
    qm_plot_signal = qm_plot_occ - initial
    fig, axes = plt.subplots(2, 5, figsize=(17, 8.5), sharex=True)
    colors = plt.get_cmap("tab10").colors
    for orbital, axis in enumerate(axes.flat):
        axis.plot(qm_times, qm_plot_signal[:, orbital], color="black", linewidth=2.0,
                  label="QM")
        for run_index, (label, times, signal) in enumerate(plot_data):
            axis.plot(times, signal[:, orbital], linewidth=1.25,
                      color=colors[run_index % len(colors)], label=label)
        axis.axhline(0.0, color="#6b7280", linewidth=0.7, linestyle="--")
        axis.set_title(f"Orbital {orbital}")
        axis.grid(alpha=0.25)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        if orbital >= 5:
            axis.set_xlabel("time (a.u.)")
        if orbital in (0, 5):
            axis.set_ylabel("occupation change")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=min(len(labels), 4),
               frameon=False)
    fig.suptitle(
        f"Short-time SepMB accuracy against same-parameter QM, t <= {args.tmax:g} a.u.",
        y=0.99,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(args.out_dir / "short_qm_accuracy_all_orbitals.png", dpi=220,
                facecolor="white")
    plt.close(fig)

    for summary in summaries:
        print(
            f"{summary['label']}: Q_qm={summary['Q_qm']:.6g} "
            f"amp_ratio={summary['amplitude_ratio']:.6g} "
            f"cosine={summary['cosine_similarity']:.6g} "
            f"error_rms={summary['error_rms']:.6e}"
        )


if __name__ == "__main__":
    main()
