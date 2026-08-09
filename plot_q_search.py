#!/usr/bin/env python3
"""Plot strict-Q histories and per-orbital repeats for a ranked search pool."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def q_windows(row: dict[str, str]) -> tuple[np.ndarray, np.ndarray]:
    windows = []
    for key, value in row.items():
        if not key.startswith("min_active_Q_") or not value:
            continue
        start, end = map(float, key.removeprefix("min_active_Q_").split("_"))
        windows.append(((start + end) / 2.0, float(value)))
    windows.sort()
    return np.array([item[0] for item in windows]), np.array([item[1] for item in windows])


def plot_q(ranking: Path, output: Path, top: int, threshold: float) -> None:
    rows = load_rows(ranking)[:top]
    fig, axis = plt.subplots(figsize=(11.5, 6.5))
    colors = plt.get_cmap("tab10").colors
    for index, row in enumerate(rows):
        time, q = q_windows(row)
        axis.plot(
            time, q, marker="o", markersize=4.2, linewidth=1.7,
            color=colors[index % len(colors)],
            label=(f"#{row['rank']} {row['case_id']} | "
                   f"min Q={float(row['required_Q']):.2f}"),
        )
    axis.axhline(threshold, color="#B91C1C", linewidth=1.4, linestyle="--",
                 label=f"target Q > {threshold:g}")
    axis.axvline(130, color="#475569", linewidth=1.1, linestyle=":",
                 label="required through 130")
    axis.set(xlabel="time-window midpoint (a.u.)", ylabel="minimum active-orbital Q",
             title="Strict convergence Q by 10 a.u. window")
    axis.set_xlim(0, 150)
    axis.set_yscale("log")
    axis.set_ylim(bottom=0.8)
    axis.grid(True, color="#CBD5E1", linewidth=0.6, alpha=0.75)
    axis.legend(loc="upper right", fontsize=8.5, frameon=True)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


def plot_orbitals(
    ranking: Path, manifest: Path, output: Path, ntraj: int,
) -> None:
    best = load_rows(ranking)[0]
    case_id = best["case_id"]
    records = [row for row in load_rows(manifest) if row["case_id"] == case_id]
    if len(records) < 3:
        raise ValueError(f"{case_id} has only {len(records)} manifest records")

    runs = []
    time = None
    for record in records:
        path = (manifest.parent / record["data_file"]).resolve()
        data = np.loadtxt(path, comments="#")
        current_time = data[:, 0]
        if time is None:
            time = current_time
        elif not np.allclose(time, current_time):
            raise ValueError(f"time grid differs in {path}")
        runs.append(data[:, 4:14])

    occupations = np.stack(runs)
    delta = occupations - occupations[:, :1, :]
    mean = delta.mean(axis=0)
    std = delta.std(axis=0, ddof=1)
    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(2, 5, figsize=(16, 8.5), sharex=True)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.09, top=0.80,
                        hspace=0.40, wspace=0.30)
    for orbital, axis in enumerate(axes.flat):
        lower = mean[:, orbital] - std[:, orbital]
        upper = mean[:, orbital] + std[:, orbital]
        axis.fill_between(time, lower, upper, color=colors[orbital], alpha=0.18,
                          linewidth=0)
        axis.plot(time, mean[:, orbital], color=colors[orbital], linewidth=1.6)
        axis.axhline(0, color="#64748B", linewidth=0.7, linestyle="--")
        low, high = float(lower.min()), float(upper.max())
        span = max(high - low, 1e-12)
        axis.set_ylim(low - 0.08 * span, high + 0.08 * span)
        axis.set_xlim(float(time[0]), float(time[-1]))
        axis.set_title(
            f"Orbital {orbital} | peak |delta n|={np.abs(mean[:, orbital]).max():.2e}",
            fontsize=9.5,
        )
        axis.grid(True, color="#CBD5E1", linewidth=0.55, alpha=0.65)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 2))
        axis.yaxis.set_major_formatter(formatter)
        if orbital >= 5:
            axis.set_xlabel("time (a.u.)")
        if orbital in (0, 5):
            axis.set_ylabel("occupation change")

    fig.suptitle(
        f"Best search result: wc={float(best['wc_eV']):g} eV, "
        f"eta={float(best['eta']):.2e}", fontsize=16, y=0.97,
    )
    fig.text(
        0.5, 0.925,
        f"{ntraj:,} trajectories x {len(runs)} repeats | solid: mean | band: +/- 1 SD | "
        f"strict Q through 130={float(best['required_Q']):.2f}",
        ha="center", fontsize=10.5, color="#334155",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top", type=int, default=6)
    parser.add_argument("--threshold", type=float, default=4.0)
    parser.add_argument("--ntraj", type=int, default=8_000_000)
    args = parser.parse_args()
    plot_q(args.ranking, args.output_dir / "top_q_by_window.png", args.top, args.threshold)
    plot_orbitals(
        args.ranking, args.manifest, args.output_dir / "best_all_orbitals.png", args.ntraj,
    )
    print(args.output_dir)


if __name__ == "__main__":
    main()
