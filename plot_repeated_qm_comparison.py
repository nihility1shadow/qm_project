#!/usr/bin/env python3
"""Plot the mean and repeat spread of SepMB runs against QM."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_short_qm_accuracy import interpolate_occupations, load_dat


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--run", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--norb", type=int, default=10)
    parser.add_argument("--nel", type=int, default=5)
    parser.add_argument("--tmax", type=float, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--subtitle", required=True)
    args = parser.parse_args()

    qm = load_dat(args.qm)
    qm = qm[qm[:, 0] <= args.tmax + 1.0e-12]
    runs = [load_dat(path) for path in args.run]
    runs = [run[run[:, 0] <= args.tmax + 1.0e-12] for run in runs]
    times = runs[0][:, 0]
    for run in runs[1:]:
        if not np.allclose(run[:, 0], times):
            raise ValueError("repeat time grids do not match")

    initial = qm[0, 4 : 4 + args.norb]
    qm_occ = interpolate_occupations(qm, times, args.norb)
    qm_signal = qm_occ - initial
    repeated_signal = np.stack(
        [run[:, 4 : 4 + args.norb] - initial for run in runs], axis=0
    )
    mean = np.mean(repeated_signal, axis=0)
    spread = np.std(repeated_signal, axis=0, ddof=1) if len(runs) > 1 else np.zeros_like(mean)

    fig, axes = plt.subplots(2, 5, figsize=(17, 8.8), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        axis.fill_between(
            times,
            mean[:, orbital] - spread[:, orbital],
            mean[:, orbital] + spread[:, orbital],
            color="#c44e52",
            alpha=0.18,
            linewidth=0,
            label="repeat SD" if orbital == 0 else None,
        )
        axis.plot(times, qm_signal[:, orbital], color="black", linewidth=2.0,
                  label="QM")
        axis.plot(times, mean[:, orbital], color="#c44e52", linewidth=1.5,
                  label="SepMB mean")
        axis.axhline(0.0, color="#6b7280", linewidth=0.7, linestyle="--")
        axis.set_title(f"Orbital {orbital}")
        axis.grid(alpha=0.25)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        if orbital >= 5:
            axis.set_xlabel("time (a.u.)")
        if orbital in (0, 5):
            axis.set_ylabel("occupation change")

    handles, labels = axes.flat[0].get_legend_handles_labels()
    order = [labels.index(name) for name in ("QM", "SepMB mean", "repeat SD")]
    fig.legend(
        [handles[index] for index in order],
        [labels[index] for index in order],
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.928),
    )
    fig.suptitle(args.title, y=0.995, fontsize=15)
    fig.text(0.5, 0.955, args.subtitle, ha="center", fontsize=10.5, color="#374151")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.89))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, facecolor="white")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
