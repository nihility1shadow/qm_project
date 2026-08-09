#!/usr/bin/env python3
"""Plot a wc-eta-Q surface from a strict convergence ranking CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_grid(path: Path, metric: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"ranking CSV is empty: {path}")
    if metric not in rows[0]:
        raise ValueError(f"metric not found in ranking CSV: {metric}")

    wc_values = np.array(sorted({float(row["wc_eV"]) for row in rows}))
    eta_values = np.array(sorted({float(row["eta"]) for row in rows}))
    q = np.full((len(eta_values), len(wc_values)), np.nan)
    wc_index = {value: index for index, value in enumerate(wc_values)}
    eta_index = {value: index for index, value in enumerate(eta_values)}
    for row in rows:
        wc = float(row["wc_eV"])
        eta = float(row["eta"])
        q[eta_index[eta], wc_index[wc]] = float(row[metric])
    return wc_values, eta_values, q


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metric", default="min_active_Q_40_50")
    parser.add_argument("--target-q", type=float, default=4.0)
    args = parser.parse_args()

    wc, eta, q = read_grid(args.ranking.resolve(), args.metric)
    x, y = np.meshgrid(wc, np.log10(eta))
    maximum = np.nanargmax(q)
    eta_index, wc_index = np.unravel_index(maximum, q.shape)
    max_wc = wc[wc_index]
    max_eta = eta[eta_index]
    max_q = q[eta_index, wc_index]

    fig = plt.figure(figsize=(16, 7.2))
    fig.subplots_adjust(
        left=0.05, right=0.96, bottom=0.14, top=0.90, wspace=0.28
    )
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    surface = ax3d.plot_surface(
        x, y, q, cmap="viridis", edgecolor="#334155", linewidth=0.35,
        antialiased=True, alpha=0.94,
    )
    ax3d.scatter(
        [max_wc], [np.log10(max_eta)], [max_q], color="#dc2626", s=55,
        depthshade=False,
    )
    ax3d.set_xlabel("wc (eV)")
    ax3d.set_ylabel("eta")
    ax3d.set_zlabel("minimum active-orbital Q (40-50 a.u.)")
    ax3d.set_yticks(np.log10(eta), [f"{value:.1e}" for value in eta])
    ax3d.tick_params(axis="y", labelsize=8)
    ax3d.view_init(elev=27, azim=-132)
    ax3d.set_title(
        f"Grid surface; max Q={max_q:.3f} at wc={max_wc:g} eV, eta={max_eta:.1e}"
    )
    fig.colorbar(surface, ax=ax3d, shrink=0.72, pad=0.08, label="Q(40-50)")

    ax2d = fig.add_subplot(1, 2, 2)
    image = ax2d.imshow(q, origin="lower", aspect="auto", cmap="viridis")
    ax2d.set_xticks(np.arange(len(wc)), [f"{value:g}" for value in wc])
    ax2d.set_yticks(np.arange(len(eta)), [f"{value:.1e}" for value in eta])
    ax2d.set_xlabel("wc (eV)")
    ax2d.set_ylabel("eta")
    ax2d.set_title(f"Top view; target Q > {args.target_q:g}")
    for row in range(q.shape[0]):
        for column in range(q.shape[1]):
            value = q[row, column]
            color = "white" if value < 0.55 * np.nanmax(q) else "black"
            ax2d.text(column, row, f"{value:.2f}", ha="center", va="center",
                      fontsize=7, color=color)
    ax2d.scatter([wc_index], [eta_index], marker="s", facecolors="none",
                 edgecolors="#dc2626", linewidths=2, s=260)
    fig.colorbar(image, ax=ax2d, shrink=0.86, label="Q(40-50)")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=210)
    plt.close(fig)
    print(
        f"output={args.output} max_wc={max_wc:g} max_eta={max_eta:.8g} "
        f"max_q={max_q:.8g}"
    )


if __name__ == "__main__":
    main()
