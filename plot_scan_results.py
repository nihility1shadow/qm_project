#!/usr/bin/env python3
"""Plot ranked convergence metrics and the best trajectory ensemble."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


COLORS = [
    "#2563eb", "#0f766e", "#b45309", "#7c3aed", "#be123c",
    "#0891b2", "#4d7c0f", "#c2410c", "#475569", "#a21caf",
]


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def case_files(manifest: Path, case_id: str) -> list[Path]:
    unique: dict[str, Path] = {}
    for row in read_rows(manifest):
        if row["case_id"] != case_id:
            continue
        path = (manifest.parent / row["data_file"]).resolve()
        unique.setdefault(file_hash(path), path)
    if not unique:
        raise ValueError(f"case not found in manifest: {case_id}")
    return list(unique.values())


def plot_ranking(rows: list[dict[str, str]], output: Path, top: int) -> None:
    selected = rows[: min(top, len(rows))]
    q_fields = [
        name for name in rows[0] if name.startswith("min_active_Q_")
    ]
    if not q_fields:
        q_fields = [
            name for name in rows[0]
            if name.startswith("Q_") and not name.startswith("Q_orb")
        ]
        prefix = "Q_"
    else:
        prefix = "min_active_Q_"
    labels = [name.removeprefix(prefix).replace("_", "-") for name in q_fields]
    x = np.arange(len(q_fields))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
    for row in selected:
        axes[0].plot(
            x,
            [float(row[name]) for name in q_fields],
            marker="o",
            linewidth=1.7,
            label=row["case_id"],
        )
    axes[0].axhline(2.0, color="#b91c1c", linestyle="--", linewidth=1)
    axes[0].set_xticks(x, labels, rotation=35, ha="right")
    axes[0].set_ylabel("minimum active-orbital Q")
    axes[0].set_title("Strict per-orbital visibility")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8)

    scatter = axes[1].scatter(
        [float(row["wc_eV"]) for row in rows],
        [float(row["eta"]) for row in rows],
        c=[float(row["required_Q"]) for row in rows],
        s=[110 if row["passes_required_interval"] == "True" else 55 for row in rows],
        cmap="viridis",
        edgecolor="#111827",
        linewidth=0.5,
    )
    for row in selected:
        axes[1].annotate(
            str(row["rank"]),
            (float(row["wc_eV"]), float(row["eta"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    axes[1].set_xlabel("wc (eV)")
    axes[1].set_ylabel("eta")
    axes[1].set_title("Required-interval Q")
    axes[1].grid(alpha=0.2)
    fig.colorbar(scatter, ax=axes[1], label="minimum Q")
    fig.savefig(output, dpi=190)
    plt.close(fig)


def plot_orbitals(
    files: list[Path], case_id: str, required_until: float, output: Path
) -> None:
    arrays = [np.loadtxt(path, comments="#") for path in files]
    time = arrays[0][:, 0]
    occupations = np.stack([array[:, 4:] for array in arrays])
    delta = occupations - occupations[:, :1, :]
    mean = delta.mean(axis=0)
    if len(arrays) >= 2:
        spread = delta.std(axis=0, ddof=1)
    else:
        spread = np.zeros_like(mean)

    norb = mean.shape[1]
    ncols = min(5, norb)
    nrows = math.ceil(norb / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 4.4 * nrows), sharex=True,
        constrained_layout=True, squeeze=False,
    )
    for orbital, ax in enumerate(axes.flat[:norb]):
        for run in delta:
            ax.plot(time, run[:, orbital], color="#94a3b8", alpha=0.34, linewidth=0.9)
        ax.fill_between(
            time,
            mean[:, orbital] - spread[:, orbital],
            mean[:, orbital] + spread[:, orbital],
            color=COLORS[orbital % len(COLORS)],
            alpha=0.16,
            linewidth=0,
        )
        ax.plot(
            time, mean[:, orbital], color=COLORS[orbital % len(COLORS)], linewidth=2.1
        )
        ax.axhline(0.0, color="#334155", linewidth=0.7, alpha=0.7)
        ax.axvline(required_until, color="#b91c1c", linestyle="--", linewidth=0.9)
        ax.set_title(f"Orbital {orbital}")
        ax.grid(alpha=0.18)
        if orbital // ncols == nrows - 1:
            ax.set_xlabel("time (a.u.)")
        if orbital % ncols == 0:
            ax.set_ylabel(r"$n_j(t)-n_j(0)$")
    for ax in axes.flat[norb:]:
        ax.set_visible(False)
    fig.suptitle(
        f"{case_id}; gray={len(arrays)} independent runs; color=mean +/- 1 SD",
        fontsize=15,
    )
    fig.savefig(output, dpi=190)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--ranking", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--top", type=int, default=6)
    args = parser.parse_args()

    rows = read_rows(args.ranking.resolve())
    if not rows:
        raise ValueError("ranking file is empty")
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    plot_ranking(rows, args.out_prefix.with_name(args.out_prefix.name + "_ranking.png"), args.top)

    best = rows[0]
    files = case_files(args.manifest.resolve(), best["case_id"])
    plot_orbitals(
        files,
        best["case_id"],
        float(best["required_until"]),
        args.out_prefix.with_name(args.out_prefix.name + "_best_orbitals.png"),
    )
    print(f"ranking_plot={args.out_prefix.with_name(args.out_prefix.name + '_ranking.png')}")
    print(f"orbital_plot={args.out_prefix.with_name(args.out_prefix.name + '_best_orbitals.png')}")


if __name__ == "__main__":
    main()
