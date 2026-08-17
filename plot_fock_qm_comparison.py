#!/usr/bin/env python3
"""Plot all-orbital and active-orbital errors for Fock/QM comparisons."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--fock", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    qm = np.loadtxt(args.qm)
    fock = np.loadtxt(args.fock)
    if qm.shape != fock.shape or not np.allclose(qm[:, 0], fock[:, 0]):
        raise ValueError("QM and Fock files must use the same time grid")

    times = qm[:, 0]
    qm_occupations = qm[:, 4:14]
    fock_occupations = fock[:, 4:14]
    qm_change = qm_occupations - qm_occupations[0]
    fock_change = fock_occupations - qm_occupations[0]
    orbital_metrics = pd.read_csv(args.metrics)
    minimum_q = orbital_metrics.groupby("orbital")["Q"].min()

    figure, axes = plt.subplots(5, 2, figsize=(13, 15), sharex=True)
    for orbital, axis in enumerate(axes.ravel()):
        axis.plot(times, qm_change[:, orbital], color="#202124", lw=1.25,
                  label="QM")
        axis.plot(times, fock_change[:, orbital], color="#d1495b", lw=0.9,
                  ls="--", label="Fock resummation")
        span = max(np.ptp(qm_change[:, orbital]),
                   np.ptp(fock_change[:, orbital]), 1.0e-12)
        lower = min(qm_change[:, orbital].min(),
                    fock_change[:, orbital].min()) - 0.08 * span
        upper = max(qm_change[:, orbital].max(),
                    fock_change[:, orbital].max()) + 0.08 * span
        axis.set_ylim(lower, upper)
        q_label = minimum_q.get(orbital)
        title_suffix = (f"min Q={q_label:.2f}" if q_label is not None
                        else "inactive, absolute scale")
        axis.set_title(f"Orbital {orbital}   {title_suffix}", fontsize=10)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        axis.grid(alpha=0.2)
        axis.set_ylabel(r"$n_j(t)-n_j(0)$")
    for axis in axes[-1]:
        axis.set_xlabel("Time (a.u.)")
    axes[0, 0].legend(loc="best", fontsize=8)
    figure.suptitle(
        "4000 a.u.: deterministic Fock resummation vs grid QM", fontsize=14
    )
    figure.tight_layout(rect=(0, 0, 1, 0.98))
    figure.savefig(args.out_dir / "fock_vs_qm_all_orbitals.png", dpi=180)
    plt.close(figure)

    active = [0, 5, 6, 7, 8, 9]
    figure, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    for orbital, axis in zip(active, axes.ravel()):
        axis.plot(times, fock_occupations[:, orbital] - qm_occupations[:, orbital],
                  color="#2a6f97", lw=0.8)
        axis.axhline(0.0, color="black", lw=0.5)
        axis.set_title(f"Orbital {orbital} absolute error")
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        axis.grid(alpha=0.2)
        axis.set_ylabel("Fock - QM")
    for axis in axes[-1]:
        axis.set_xlabel("Time (a.u.)")
    figure.tight_layout()
    figure.savefig(args.out_dir / "fock_active_orbital_error.png", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
