#!/usr/bin/env python3
"""Estimate SepMB trajectory requirements from repeated long-time runs."""

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


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(signal: np.ndarray, error: np.ndarray) -> float:
    error_rms = rms(error)
    return rms(signal) / error_rms if error_rms > 0.0 else float("inf")


def parse_run(value: str) -> dict[str, object]:
    fields = value.split("=", 1)
    if len(fields) != 2:
        raise argparse.ArgumentTypeError(
            "run must use LABEL,NTRAJ,BACK_REPLICAS,SECONDS=PATH"
        )
    metadata, raw_path = fields
    parts = metadata.split(",")
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            "run metadata must contain LABEL,NTRAJ,BACK_REPLICAS,SECONDS"
        )
    return {
        "label": parts[0],
        "ntraj": int(parts[1]),
        "back_replicas": int(parts[2]),
        "seconds": float(parts[3]),
        "path": Path(raw_path),
    }


def parse_orbitals(value: str) -> list[int]:
    return [int(item) for item in value.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--run", required=True, action="append", type=parse_run)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--target-q", type=float, default=10.0)
    parser.add_argument("--norb", type=int, default=10)
    parser.add_argument("--active-orbitals", type=parse_orbitals, default=[0, 5, 6, 7, 8, 9])
    parser.add_argument("--window-edges", type=float, nargs="+", default=[0, 50, 100, 150, 200])
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    qm = load_dat(args.qm)
    initial = qm[0, 4 : 4 + args.norb]
    active = np.asarray(args.active_orbitals, dtype=int)
    rows: list[dict[str, object]] = []
    loaded_runs: list[tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray]] = []

    for spec in args.run:
        data = load_dat(spec["path"])
        times = data[:, 0]
        qm_occ = interpolate_occupations(qm, times, args.norb)
        signal = qm_occ - initial
        simulation = data[:, 4 : 4 + args.norb] - initial
        error = simulation - signal
        loaded_runs.append((spec, times, signal, simulation))

        for start, stop in zip(args.window_edges[:-1], args.window_edges[1:]):
            mask = (times >= start - 1.0e-12) & (times <= stop + 1.0e-12)
            for scope, columns in (("all", np.arange(args.norb)), ("active", active)):
                q = quality(signal[mask][:, columns], error[mask][:, columns])
                multiplier = (args.target_q / q) ** 2 if q > 0.0 else float("inf")
                rows.append(
                    {
                        "label": spec["label"],
                        "scope": scope,
                        "window_start": start,
                        "window_stop": stop,
                        "ntraj": spec["ntraj"],
                        "back_replicas": spec["back_replicas"],
                        "runtime_seconds_64_cores": spec["seconds"],
                        "Q": q,
                        "qm_signal_rms": rms(signal[mask][:, columns]),
                        "error_rms": rms(error[mask][:, columns]),
                        "target_Q": args.target_q,
                        "trajectory_multiplier": multiplier,
                        "estimated_ntraj_for_target": spec["ntraj"] * multiplier,
                        "estimated_runtime_seconds_64_cores": spec["seconds"] * multiplier,
                    }
                )

    fieldnames = list(rows[0].keys())
    with (args.out_dir / "t200_sampling_requirement.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    (args.out_dir / "t200_sampling_requirement.json").write_text(
        json.dumps(rows, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    # Show how each fixed-budget allocation behaves in consecutive windows.
    active_rows = [row for row in rows if row["scope"] == "active"]
    labels = list(dict.fromkeys(str(row["label"]) for row in active_rows))
    windows = [
        f"{start:g}-{stop:g}"
        for start, stop in zip(args.window_edges[:-1], args.window_edges[1:])
    ]
    x = np.arange(len(windows))
    fig, axis = plt.subplots(figsize=(10.5, 5.8))
    for label in labels:
        selected = [row for row in active_rows if row["label"] == label]
        axis.plot(x, [row["Q"] for row in selected], marker="o", linewidth=1.6, label=label)
    axis.axhline(args.target_q, color="black", linestyle="--", linewidth=1.0,
                 label=f"target Q={args.target_q:g}")
    axis.set_xticks(x, windows)
    axis.set_xlabel("time window (a.u.)")
    axis.set_ylabel("Q against same-parameter QM (active orbitals)")
    axis.set_yscale("log")
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(args.out_dir / "t200_window_q.png", dpi=220, facecolor="white")
    plt.close(fig)

    for row in active_rows:
        print(
            f"{row['label']} {row['window_start']:g}-{row['window_stop']:g}: "
            f"Q={row['Q']:.4g}, N(Q={args.target_q:g})="
            f"{row['estimated_ntraj_for_target']:.4g}, "
            f"runtime64={row['estimated_runtime_seconds_64_cores'] / 3600.0:.3g} h"
        )


if __name__ == "__main__":
    main()
