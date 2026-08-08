#!/usr/bin/env python3
"""Rank O2-HOMO Poisson scans and plot all orbital occupations.

This script is intentionally local-only.  It reads one scan directory whose
children each contain output.txt and ahm-sepmb-s10-n5-*.dat.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PARAM_RE = re.compile(
    r"wc_eV=(?P<wc>[-+0-9.eE]+).*eta=(?P<eta>[-+0-9.eE]+)"
    r".*delE_eV=(?P<dele>[-+0-9.eE]+).*nstep=(?P<nstep>\d+)"
)


def moving_average(values: np.ndarray, window: int = 7) -> np.ndarray:
    if values.shape[0] < window:
        return values.copy()
    pad = window // 2
    padded = np.pad(values, ((pad, pad), (0, 0)), mode="reflect")
    kernel = np.ones(window, dtype=float) / window
    return np.column_stack(
        [np.convolve(padded[:, j], kernel, mode="valid") for j in range(values.shape[1])]
    )


def valid_until(time: np.ndarray, valid: np.ndarray) -> float:
    failed = np.flatnonzero(~valid)
    if not failed.size:
        return float(time[-1])
    if failed[0] == 0:
        return 0.0
    return float(time[failed[0] - 1])


def read_candidate(directory: Path, noise_limit: float) -> dict:
    output_path = directory / "output.txt"
    if not output_path.exists():
        output_path = directory / "program.out"
    output = output_path.read_text(encoding="utf-8", errors="replace")
    match = PARAM_RE.search(output)
    if not match:
        raise ValueError(f"missing #AHAU_PARAMS in {output_path}")

    data_files = sorted(directory.glob("ahm-sepmb-s*-n*-*.dat"))
    if len(data_files) != 1:
        raise ValueError(f"expected one SepMB data file in {directory}, found {len(data_files)}")
    raw = np.loadtxt(data_files[0], comments="#")
    ntraj_match = re.search(r"-(\d+)\.dat$", data_files[0].name)
    if not ntraj_match:
        raise ValueError(f"cannot read ntraj from {data_files[0].name}")
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    if raw.shape[1] < 5:
        raise ValueError(f"not enough columns in {data_files[0]}")

    time = raw[:, 0]
    nel = raw[:, 1]
    occ = raw[:, 4:]
    smooth = moving_average(occ)
    initial = occ[0]
    smooth_delta = smooth - initial
    raw_delta = occ - initial
    trace_error = occ.sum(axis=1) - nel

    trend_t = np.sqrt(np.mean(smooth_delta * smooth_delta, axis=1))
    raw_t = np.sqrt(np.mean(raw_delta * raw_delta, axis=1))
    noise = occ - smooth
    noise_t = np.sqrt(np.mean(noise * noise, axis=1))

    physical = (
        (occ.min(axis=1) >= -0.25)
        & (occ.max(axis=1) <= 1.25)
        & (np.abs(trace_error) <= 0.25)
    )
    stable_t = valid_until(time, physical)
    precision_t = valid_until(time, noise_t <= noise_limit)
    precision_mask = time <= precision_t

    trend_rms = float(np.sqrt(np.mean(smooth_delta * smooth_delta)))
    trend_peak = float(np.max(trend_t))
    raw_peak = float(np.max(raw_t))
    noise_rms = float(np.sqrt(np.mean(noise * noise)))
    trace_rmse = float(np.sqrt(np.mean(trace_error * trace_error)))
    trace_max = float(np.max(np.abs(trace_error)))
    roughness = float(np.sqrt(np.mean(np.diff(occ, n=2, axis=0) ** 2)))
    physical_fraction = float(np.mean(physical))

    trend_before_precision = float(np.max(trend_t[precision_mask]))

    return {
        "name": directory.name,
        "directory": directory,
        "data_file": data_files[0],
        "wc_eV": float(match.group("wc")),
        "eta": float(match.group("eta")),
        "delE_eV": float(match.group("dele")),
        "nstep": int(match.group("nstep")),
        "ntraj": int(ntraj_match.group(1)),
        "n_time": int(time.size),
        "tmax": float(time[-1]),
        "stable_t": stable_t,
        "precision_t": precision_t,
        "noise_limit": noise_limit,
        "physical_fraction": physical_fraction,
        "occ_min": float(occ.min()),
        "occ_max": float(occ.max()),
        "trace_rmse": trace_rmse,
        "trace_max": trace_max,
        "trend_rms": trend_rms,
        "trend_peak": trend_peak,
        "trend_before_precision": trend_before_precision,
        "raw_peak": raw_peak,
        "noise_rms": noise_rms,
        "roughness": roughness,
        "time": time,
        "occ": occ,
        "trend_t": trend_t,
        "noise_t": noise_t,
        "trace_error": trace_error,
    }


def write_metrics(candidates: list[dict], root: Path) -> None:
    fields = [
        "rank", "name", "wc_eV", "eta", "delE_eV", "nstep", "ntraj", "n_time", "tmax",
        "stable_t", "precision_t", "noise_limit", "physical_fraction", "occ_min",
        "occ_max", "trace_rmse", "trace_max", "trend_rms", "trend_peak",
        "trend_before_precision", "raw_peak", "noise_rms", "roughness",
    ]
    rows = []
    for rank, candidate in enumerate(candidates, 1):
        row = {field: candidate.get(field) for field in fields}
        row["rank"] = rank
        rows.append(row)

    with (root / "scan_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (root / "scan_metrics.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def plot_overview(candidates: list[dict], root: Path) -> None:
    tmax_label = f"{max(c['tmax'] for c in candidates):g}".replace(".", "p")
    fig, ax = plt.subplots(figsize=(11, 6.2), constrained_layout=True)
    for rank, candidate in enumerate(candidates, 1):
        ax.plot(
            candidate["time"], candidate["trend_t"], linewidth=1.3,
            label=f"{rank}. {candidate['name']}",
        )
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"orbital RMS change $\sqrt{\langle[n_j(t)-n_j(0)]^2\rangle_j}$")
    ax.set_xlim(0, max(c["tmax"] for c in candidates))
    ax.grid(alpha=0.22)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.savefig(root / f"trend-overview-0-{tmax_label}.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 6.3), constrained_layout=True)
    precision_times = np.array([candidate["precision_t"] for candidate in candidates])
    max_visible = max(candidate["trend_before_precision"] for candidate in candidates)
    sizes = 90.0 + 600.0 * np.array(
        [candidate["trend_before_precision"] / max_visible for candidate in candidates]
    )
    scatter = ax.scatter(
        [candidate["wc_eV"] for candidate in candidates],
        [candidate["eta"] for candidate in candidates],
        c=precision_times, s=sizes, cmap="viridis", edgecolors="black", linewidths=0.5,
    )
    for rank, candidate in enumerate(candidates, 1):
        ax.annotate(str(rank), (candidate["wc_eV"], candidate["eta"]),
                    xytext=(5, 5), textcoords="offset points", fontsize=9)
    ax.set_xlabel(r"$\omega_c$ (eV)")
    ax.set_ylabel(r"$\eta$")
    ax.set_yscale("log")
    ax.grid(alpha=0.2)
    colorbar = fig.colorbar(scatter, ax=ax)
    colorbar.set_label("converged window (a.u.)")
    fig.savefig(root / f"parameter-map-0-{tmax_label}.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.6), constrained_layout=True)
    for rank, candidate in enumerate(candidates, 1):
        ax.plot(candidate["time"], candidate["trace_error"], linewidth=1.1,
                label=f"{rank}. {candidate['name']}")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$\sum_j n_j-N_{el}$")
    ax.grid(alpha=0.22)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.savefig(root / f"electron-number-error-0-{tmax_label}.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.6), constrained_layout=True)
    for rank, candidate in enumerate(candidates, 1):
        ax.plot(candidate["time"], candidate["noise_t"], linewidth=1.1,
                label=f"{rank}. {candidate['name']}")
    ax.axhline(candidates[0]["noise_limit"], color="black", linewidth=0.9,
               linestyle="--", label="noise limit")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("local high-frequency occupation noise")
    ax.grid(alpha=0.22)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.savefig(root / f"noise-overview-0-{tmax_label}.png", dpi=180)
    plt.close(fig)


def plot_orbitals(candidate: dict, root: Path) -> None:
    time = candidate["time"]
    occ = candidate["occ"]
    n_orb = occ.shape[1]
    ncols = 5
    nrows = math.ceil(n_orb / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(15, 3.7 * nrows), sharex=True, constrained_layout=True
    )
    axes = np.atleast_1d(axes).ravel()
    for j, ax in enumerate(axes):
        if j >= n_orb:
            ax.set_visible(False)
            continue
        y = occ[:, j]
        baseline = y[0]
        lo, hi = np.quantile(y, [0.01, 0.99])
        span = max(float(hi - lo), 2.0e-3)
        margin = 0.15 * span
        ax.plot(time, y, color="tab:blue", linewidth=1.15)
        ax.axhline(baseline, color="black", linewidth=0.7, alpha=0.5)
        ax.set_ylim(float(lo - margin), float(hi + margin))
        ax.set_title(f"orbital {j}")
        ax.grid(alpha=0.2)
        if j % ncols == 0:
            ax.set_ylabel(r"$n_j(t)$")
        if j // ncols == nrows - 1:
            ax.set_xlabel("time (a.u.)")
    fig.suptitle(
        f"{candidate['name']}: wc={candidate['wc_eV']:.4g} eV, "
        f"eta={candidate['eta']:.3g}, ntraj={candidate['ntraj']}",
        fontsize=13,
    )
    fig.savefig(root / f"all-orbitals-{candidate['name']}.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scan_dir", type=Path)
    parser.add_argument("--top", type=int, default=4)
    parser.add_argument("--noise-limit", type=float, default=2.0e-3)
    args = parser.parse_args()
    root = args.scan_dir.resolve()

    candidates = [
        read_candidate(path, args.noise_limit)
        for path in sorted(root.iterdir())
        if path.is_dir()
        and ((path / "output.txt").exists() or (path / "program.out").exists())
    ]
    if not candidates:
        raise SystemExit(f"no candidates found under {root}")
    candidates.sort(
        key=lambda candidate: (
            candidate["precision_t"], candidate["trend_before_precision"],
            -candidate["noise_rms"],
        ),
        reverse=True,
    )

    write_metrics(candidates, root)
    plot_overview(candidates, root)
    for candidate in candidates[: args.top]:
        plot_orbitals(candidate, root)

    print(
        "rank,name,wc_eV,eta,stable_t,precision_t,trend_before_precision,"
        "trend_peak,noise_rms,trace_rmse"
    )
    for rank, candidate in enumerate(candidates, 1):
        print(
            f"{rank},{candidate['name']},{candidate['wc_eV']:.6g},"
            f"{candidate['eta']:.8g},{candidate['stable_t']:.3f},"
            f"{candidate['precision_t']:.3f},{candidate['trend_before_precision']:.6g},"
            f"{candidate['trend_peak']:.6g},{candidate['noise_rms']:.6g},"
            f"{candidate['trace_rmse']:.6g}"
        )


if __name__ == "__main__":
    main()
