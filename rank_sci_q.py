#!/usr/bin/env python3
"""Combine SCI/TCI output with independent-repeat signal/noise checks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SEGMENTS = ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50))


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


def q_metrics(runs: list[np.ndarray]) -> dict:
    time = runs[0][:, 0]
    occupations = np.stack([run[:, 4:14] for run in runs])
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

    q_values = []
    for start, end in SEGMENTS:
        mask = (time >= start) & (time <= end)
        signal = np.sqrt(np.mean(mean[mask] ** 2))
        noise = np.sqrt(np.mean(variance[mask]))
        q_values.append(float(signal / max(noise, 1.0e-30)))

    return {
        "q": q_values,
        "signal_rms": float(np.sqrt(np.mean(mean**2))),
        "signal_peak": float(np.max(np.sqrt(np.mean(mean**2, axis=1)))),
        "noise_rms": float(np.sqrt(np.mean(variance))),
        "particle_number_max_error": float(
            np.max(np.abs(occupations.sum(axis=2) - occupations[:, :, 0] * 0 - 5.0))
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sci-json", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--min-sci", type=float, default=80.0)
    parser.add_argument("--min-q", type=float, default=2.0)
    args = parser.parse_args()

    cases = read_manifest(args.manifest.resolve())
    sci_document = json.loads(args.sci_json.read_text(encoding="utf-8"))
    sci_cases = {case["case_id"]: case for case in sci_document["cases"]}
    rows = []

    for case_id, case in cases.items():
        runs = independent_runs(case["files"])
        metrics = q_metrics(runs)
        sci = sci_cases[case_id]
        segments = sci["segments"]
        sci_30 = next(segment["SCI"] for segment in segments if 30 <= segment["start"] < 40)
        sci_40 = next(segment["SCI"] for segment in segments if 40 <= segment["start"])
        tail_30 = min(sci_30, sci_40)
        passes = tail_30 >= args.min_sci and metrics["q"][3] >= args.min_q
        row = {
            "case_id": case_id,
            "wc_eV": case["wc_eV"],
            "eta": case["eta"],
            "ntraj_per_run": case["ntraj"],
            "independent_repeats": len(runs),
            "status": sci["status"],
            "TCI": sci["total_index"],
            "SCI_30_40": sci_30,
            "SCI_40_50": sci_40,
            "tail_SCI_from_30": tail_30,
            "Q_0_10": metrics["q"][0],
            "Q_10_20": metrics["q"][1],
            "Q_20_30": metrics["q"][2],
            "Q_30_40": metrics["q"][3],
            "Q_40_50": metrics["q"][4],
            "signal_rms": metrics["signal_rms"],
            "signal_peak": metrics["signal_peak"],
            "repeat_noise_rms": metrics["noise_rms"],
            "particle_number_max_error": metrics["particle_number_max_error"],
            "passes_40t": passes,
        }
        rows.append(row)

    rows.sort(
        key=lambda row: (
            row["passes_40t"],
            row["tail_SCI_from_30"],
            row["TCI"],
            row["Q_30_40"],
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

    top = rows[: min(8, len(rows))]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    ax = axes[0, 0]
    scatter = ax.scatter(
        [row["wc_eV"] for row in rows],
        [row["eta"] for row in rows],
        c=[row["tail_SCI_from_30"] for row in rows],
        s=[120 if row["passes_40t"] else 55 for row in rows],
        cmap="viridis", edgecolor="black", linewidth=0.5,
    )
    ax.set_xlabel("wc (eV)")
    ax.set_ylabel("eta")
    ax.set_title("Tail SCI from 30t")
    fig.colorbar(scatter, ax=ax, label="SCI")

    x = np.arange(len(SEGMENTS))
    labels = [f"{start}-{end}" for start, end in SEGMENTS]
    ax = axes[0, 1]
    for row in top:
        sci = sci_cases[row["case_id"]]
        ax.plot(x, [segment["SCI"] for segment in sci["segments"]], marker="o", label=row["case_id"])
    ax.axhline(args.min_sci, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(x, labels)
    ax.set_ylabel("SCI")
    ax.set_title("SCI by time segment")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=7, ncol=2)

    ax = axes[1, 0]
    for row in top:
        ax.plot(x, [row[f"Q_{label.replace('-', '_')}"] for label in labels], marker="o", label=row["case_id"])
    ax.axhline(args.min_q, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(x, labels)
    ax.set_ylabel("repeat signal/noise Q")
    ax.set_title("Independent-repeat visibility")
    ax.grid(alpha=0.2)

    ax = axes[1, 1]
    for row in rows:
        color = "#16803c" if row["passes_40t"] else "#64748b"
        ax.scatter(row["signal_peak"], row["TCI"], color=color, s=55)
    ax.set_xlabel("peak orbital RMS change")
    ax.set_ylabel("TCI")
    ax.set_title("Physical signal retained vs convergence")
    ax.grid(alpha=0.2)
    fig.suptitle("SCI/TCI and independent-repeat screening")
    fig.savefig(args.out_prefix.with_suffix(".png"), dpi=190)
    plt.close(fig)

    print("rank,case,pass,tail_SCI_30,Q_30_40,TCI,signal_peak,repeats")
    for row in rows[:10]:
        print(
            f"{row['rank']},{row['case_id']},{row['passes_40t']},"
            f"{row['tail_SCI_from_30']:.2f},{row['Q_30_40']:.2f},"
            f"{row['TCI']:.2f},{row['signal_peak']:.3e},{row['independent_repeats']}"
        )


if __name__ == "__main__":
    main()
