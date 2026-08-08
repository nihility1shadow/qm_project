#!/usr/bin/env python3
"""Rank stochastic trajectories by independent-repeat signal visibility."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> dict[str, dict]:
    cases: dict[str, dict] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            case = cases.setdefault(
                row["case_id"],
                {
                    "files": [],
                    "wc_eV": float(row["param.wc_eV"]),
                    "eta": float(row["param.eta"]),
                    "ntraj": int(row["param.ntraj"]),
                },
            )
            case["files"].append((path.parent / row["data_file"]).resolve())
    return cases


def independent_runs(files: list[Path]) -> list[np.ndarray]:
    unique: dict[str, Path] = {}
    for path in files:
        unique.setdefault(file_hash(path), path)
    return [np.loadtxt(path, comments="#") for path in unique.values()]


def segment_bounds(time: np.ndarray, width: float) -> list[tuple[float, float]]:
    if width <= 0:
        raise ValueError("segment width must be positive")
    start = float(time[0])
    stop = float(time[-1])
    segments = []
    while start < stop - 1.0e-12:
        end = min(start + width, stop)
        segments.append((start, end))
        start = end
    return segments


def segment_label(start: float, end: float) -> str:
    return f"{start:g}-{end:g}"


def q_metrics(runs: list[np.ndarray], segment_width: float) -> dict:
    time = runs[0][:, 0]
    occupations = np.stack([run[:, 4:] for run in runs])
    initial = occupations[:, :1, :]
    delta = occupations - initial
    mean = delta.mean(axis=0)

    if len(runs) >= 2:
        variance = delta.var(axis=0, ddof=1)
    else:
        smooth = np.column_stack(
            [np.convolve(mean[:, j], np.ones(7) / 7, mode="same") for j in range(mean.shape[1])]
        )
        variance = (mean - smooth) ** 2

    segments = segment_bounds(time, segment_width)
    q_values = {}
    for start, end in segments:
        mask = (time >= start) & (time <= end)
        signal = np.sqrt(np.mean(mean[mask] ** 2))
        noise = np.sqrt(np.mean(variance[mask]))
        q_values[segment_label(start, end)] = float(signal / max(noise, 1.0e-30))

    return {
        "q": q_values,
        "segments": segments,
        "time_end": float(time[-1]),
        "signal_rms": float(np.sqrt(np.mean(mean**2))),
        "signal_peak": float(np.max(np.sqrt(np.mean(mean**2, axis=1)))),
        "noise_rms": float(np.sqrt(np.mean(variance))),
        "particle_number_max_error": float(
            np.max(
                np.abs(
                    occupations.sum(axis=2)
                    - occupations[:, 0, :].sum(axis=1)[:, np.newaxis]
                )
            )
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument(
        "--sci-json", type=Path,
        help="optional flat-tail diagnostic; it does not decide physical convergence",
    )
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--min-q", type=float, default=2.0)
    parser.add_argument("--min-repeats", type=int, default=3)
    parser.add_argument("--segment-width", type=float, default=10.0)
    parser.add_argument(
        "--required-until", type=float,
        help="all segments up to this time must pass; defaults to the data end",
    )
    parser.add_argument(
        "--tail-width", type=float, default=20.0,
        help="final-time window used for the pass/fail decision",
    )
    args = parser.parse_args()

    cases = read_manifest(args.manifest.resolve())
    sci_cases = {}
    if args.sci_json:
        sci_document = json.loads(args.sci_json.read_text(encoding="utf-8"))
        sci_cases = {case["case_id"]: case for case in sci_document["cases"]}
    rows = []

    for case_id, case in cases.items():
        runs = independent_runs(case["files"])
        metrics = q_metrics(runs, args.segment_width)
        required_until = metrics["time_end"] if args.required_until is None else args.required_until
        if required_until > metrics["time_end"] + 1.0e-12:
            raise ValueError(
                f"{case_id}: required time {required_until:g} exceeds data end "
                f"{metrics['time_end']:g}"
            )
        required_q_labels = [
            segment_label(start, end)
            for start, end in metrics["segments"]
            if end <= required_until + 1.0e-12
        ]
        if not required_q_labels:
            raise ValueError(f"{case_id}: no complete Q segment before required time")
        required_q = min(metrics["q"][label] for label in required_q_labels)
        tail_start = max(float(metrics["segments"][0][0]), metrics["time_end"] - args.tail_width)
        tail_q_labels = [
            segment_label(start, end)
            for start, end in metrics["segments"]
            if end > tail_start
        ]
        tail_q = min(metrics["q"][label] for label in tail_q_labels)
        passes = (
            required_q >= args.min_q
            and len(runs) >= args.min_repeats
            and metrics["particle_number_max_error"] <= 1.0e-8
        )
        row = {
            "case_id": case_id,
            "wc_eV": case["wc_eV"],
            "eta": case["eta"],
            "ntraj_per_run": case["ntraj"],
            "independent_repeats": len(runs),
            "time_end": metrics["time_end"],
            "required_until": required_until,
            "required_Q": required_q,
            "tail_window_start": tail_start,
            "tail_Q": tail_q,
            "signal_rms": metrics["signal_rms"],
            "signal_peak": metrics["signal_peak"],
            "repeat_noise_rms": metrics["noise_rms"],
            "particle_number_max_error": metrics["particle_number_max_error"],
            "passes_required_interval": passes,
        }
        if case_id in sci_cases:
            sci = sci_cases[case_id]
            row["flat_tail_status"] = sci["status"]
            row["flat_tail_TCI"] = sci["total_index"]
        for label, value in metrics["q"].items():
            row[f"Q_{label.replace('-', '_')}"] = value
        rows.append(row)

    rows.sort(
        key=lambda row: (
            row["passes_required_interval"],
            row["required_Q"],
            -row["repeat_noise_rms"],
        ),
        reverse=True,
    )
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fields = ["rank"] + [key for key in rows[0] if key != "rank"]
    with args.out_prefix.with_suffix(".csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    args.out_prefix.with_suffix(".json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("rank,case,pass,required_until,required_Q,signal_peak,noise_rms,repeats")
    for row in rows[:10]:
        print(
            f"{row['rank']},{row['case_id']},{row['passes_required_interval']},"
            f"{row['required_until']:g},{row['required_Q']:.2f},"
            f"{row['signal_peak']:.3e},{row['repeat_noise_rms']:.3e},"
            f"{row['independent_repeats']}"
        )


if __name__ == "__main__":
    main()
