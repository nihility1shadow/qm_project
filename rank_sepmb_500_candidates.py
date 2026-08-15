#!/usr/bin/env python3
"""Rank long-time SepMB candidates against their matching QM references."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from analyze_short_qm_accuracy import interpolate_occupations, load_dat


DEFAULT_ACTIVE = (0, 5, 6, 7, 8, 9)


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def q_value(signal: np.ndarray, error: np.ndarray) -> float:
    error_rms = rms(error)
    return rms(signal) / error_rms if error_rms > 0.0 else math.inf


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left.ravel(), right.ravel()) / denominator) \
        if denominator > 0.0 else math.nan


def resolve_path(manifest: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else manifest.parent / path


def analyze_case(
    row: dict[str, str], manifest: Path, active: np.ndarray, norb: int,
    nel: int, tmax: float, target_stop: float, window: float, target_q: float,
) -> dict[str, object]:
    simulation_path = resolve_path(manifest, row["simulation_path"])
    qm_path = resolve_path(manifest, row["qm_path"])
    simulation = load_dat(simulation_path)
    simulation = simulation[simulation[:, 0] <= tmax + 1.0e-12]
    times = simulation[:, 0]
    qm = load_dat(qm_path)
    qm_occ = interpolate_occupations(qm, times, norb)
    initial = qm_occ[0]
    signal = qm_occ - initial
    estimate = simulation[:, 4 : 4 + norb] - initial
    error = estimate - signal

    windows: dict[str, float] = {}
    orbital_windows: dict[str, float] = {}
    start = 0.0
    while start < tmax - 1.0e-12:
        stop = min(start + window, tmax)
        mask = (times >= start - 1.0e-12) & (times <= stop + 1.0e-12)
        windows[f"Q_{start:g}_{stop:g}"] = q_value(
            signal[mask][:, active], error[mask][:, active]
        )
        orbital_windows[f"min_orbital_Q_{start:g}_{stop:g}"] = min(
            q_value(signal[mask, orbital], error[mask, orbital])
            for orbital in active
        )
        start = stop

    target_values = [
        value for key, value in windows.items()
        if float(key.split("_")[-1]) <= target_stop + 1.0e-12
    ]
    minimum_q = min(target_values)
    minimum_orbital_q = min(
        value for key, value in orbital_windows.items()
        if float(key.split("_")[-1]) <= target_stop + 1.0e-12
    )
    target_mask = times <= target_stop + 1.0e-12
    target_signal = signal[target_mask][:, active]
    target_estimate = estimate[target_mask][:, active]
    inactive = np.asarray(
        [orbital for orbital in range(norb) if orbital not in set(active)], dtype=int
    )
    ntraj = int(row["ntraj"])
    multiplier = (target_q / minimum_q) ** 2 if minimum_q > 0.0 else math.inf
    orbital_multiplier = (target_q / minimum_orbital_q) ** 2 \
        if minimum_orbital_q > 0.0 else math.inf
    runtime = float(row.get("runtime_seconds") or "nan")
    result: dict[str, object] = dict(row)
    result.update(
        {
            "n_points": int(len(times)),
            "actual_tmax": float(times[-1]),
            "min_Q_0_target": minimum_q,
            "min_active_orbital_Q_0_target": minimum_orbital_q,
            "target_stop": target_stop,
            "target_Q": target_q,
            "projected_multiplier": multiplier,
            "projected_ntraj_Q_target": ntraj * multiplier,
            "projected_ntraj_orbital_Q_target": ntraj * orbital_multiplier,
            "Q_per_sqrt_ntraj": minimum_q / math.sqrt(ntraj),
            "projected_runtime_seconds": runtime * multiplier,
            "signal_rms_active_0_target": rms(target_signal),
            "signal_peak_active_0_target": float(np.max(np.abs(target_signal))),
            "amplitude_ratio_active_0_target": (
                rms(target_estimate) / rms(target_signal)
                if rms(target_signal) > 0.0 else math.nan
            ),
            "cosine_active_0_target": cosine(target_signal, target_estimate),
            "inactive_absolute_error_rms_0_target": (
                rms(error[target_mask][:, inactive]) if len(inactive) else 0.0
            ),
            "max_particle_error": float(
                np.max(np.abs(np.sum(simulation[:, 4 : 4 + norb], axis=1) - nel))
            ),
        }
    )
    result.update(windows)
    result.update(orbital_windows)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--tmax", default=500.0, type=float)
    parser.add_argument("--target-stop", default=400.0, type=float)
    parser.add_argument("--window", default=50.0, type=float)
    parser.add_argument("--target-q", default=10.0, type=float)
    parser.add_argument("--norb", default=10, type=int)
    parser.add_argument("--nel", default=5, type=int)
    parser.add_argument("--active", default=",".join(map(str, DEFAULT_ACTIVE)))
    args = parser.parse_args()

    active = np.asarray([int(item) for item in args.active.split(",")], dtype=int)
    with args.manifest.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    results = [
        analyze_case(
            row, args.manifest, active, args.norb, args.nel, args.tmax,
            args.target_stop, args.window, args.target_q,
        )
        for row in rows
    ]
    results.sort(key=lambda item: float(item["projected_ntraj_Q_target"]))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for result in results:
        for key in result:
            if key not in fieldnames:
                fieldnames.append(key)
    with (args.out_dir / "candidate_ranking.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    (args.out_dir / "candidate_ranking.json").write_text(
        json.dumps(results, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    for rank, result in enumerate(results, 1):
        print(
            f"{rank:2d} {result['case_id']}: "
            f"minQ={float(result['min_Q_0_target']):.5g}, "
            f"projected_N={float(result['projected_ntraj_Q_target']):.5g}, "
            f"runtime={float(result['runtime_seconds']):.5g}s"
        )


if __name__ == "__main__":
    main()
