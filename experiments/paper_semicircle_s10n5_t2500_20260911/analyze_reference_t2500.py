#!/usr/bin/env python3
"""Validate bounded-reference propagation against the 2500-au grid-QM run."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
QM_PATH = ROOT / "648247" / "ahm-qm-s10-n5.dat"
CASES = {
    "D2/F384 (186 states)": ROOT / "648250" / "reference-observables.dat",
    "D3/F384 (246 states)": ROOT / "648249" / "reference-observables.dat",
    "D4/F384 (252 states)": ROOT / "648251" / "reference-observables.dat",
}
ACTIVE = [0, 5, 6, 7, 8, 9]
NORB = 10
NEL = 5
WINDOW = 100.0


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(values * values)))


def q_value(signal: np.ndarray, error: np.ndarray) -> float:
    denominator = rms(error)
    return rms(signal) / denominator if denominator > 0.0 else float("inf")


def load_reference(path: Path) -> np.ndarray:
    raw = np.loadtxt(path, comments="#")
    norm = raw[:, 1] / NEL
    result = raw.copy()
    result[:, 2:] /= norm[:, None]
    result[:, 1] = NEL
    return result


def resources(folder: Path) -> dict[str, float | int | None]:
    text = (folder / "program.out").read_text(encoding="utf-8")

    def number(pattern: str) -> float | None:
        match = re.search(pattern, text)
        return float(match.group(1)) if match else None

    state_count = number(r"active=1 states=(\d+)")
    return {
        "reference_states": int(state_count) if state_count is not None else None,
        "reference_seconds": number(r"reference_seconds=([0-9.eE+-]+)"),
        "max_rank_peak_rss_kib": number(r"max_rank_peak_rss_kb=(\d+)"),
        "mpi_ranks": int(number(r"#AHAU_MPI nproc=(\d+)") or 0),
    }


def metrics(qm: np.ndarray, estimate: np.ndarray) -> dict[str, object]:
    assert np.array_equal(qm[:, 0], estimate[:, 0])
    time = qm[:, 0]
    initial = qm[0, 4 : 4 + NORB]
    signal = qm[:, 4 : 4 + NORB] - initial
    error = estimate[:, 4 : 4 + NORB] - qm[:, 4 : 4 + NORB]
    windows: list[dict[str, object]] = []
    for start in np.arange(0.0, time[-1], WINDOW):
        stop = min(float(time[-1]), float(start + WINDOW))
        mask = (time >= start - 1.0e-12) & (time <= stop + 1.0e-12)
        per_orbital = [
            q_value(signal[mask, orbital], error[mask, orbital])
            for orbital in range(NORB)
        ]
        weakest = min(ACTIVE, key=lambda orbital: per_orbital[orbital])
        windows.append(
            {
                "start": float(start),
                "stop": stop,
                "Q_active": q_value(signal[mask][:, ACTIVE], error[mask][:, ACTIVE]),
                "Q_weakest_active": per_orbital[weakest],
                "weakest_active_orbital": weakest,
                "Q_per_orbital": per_orbital,
                "max_abs_orbital_error": float(np.max(np.abs(error[mask]))),
            }
        )
    return {
        "Q_active_global": q_value(signal[:, ACTIVE], error[:, ACTIVE]),
        "Q_min_active_window": min(row["Q_active"] for row in windows),
        "Q_min_active_orbital_window": min(
            row["Q_weakest_active"] for row in windows
        ),
        "max_abs_orbital_error": float(np.max(np.abs(error))),
        "max_particle_error": float(
            np.max(np.abs(np.sum(estimate[:, 4 : 4 + NORB], axis=1) - NEL))
        ),
        "windows": windows,
    }


def main() -> None:
    qm = np.loadtxt(QM_PATH, comments="#")
    assert len(qm) == 5001 and qm[0, 0] == 0.0 and qm[-1, 0] == 2500.0
    estimates = {label: load_reference(path) for label, path in CASES.items()}
    for estimate in estimates.values():
        assert estimate.shape == qm.shape
        assert np.isfinite(estimate).all()

    report = {
        "scope": "10 orbitals / 5 electrons / 0-2500 au / dt=0.5",
        "active_orbitals": ACTIVE,
        "Q_definition": (
            "RMS(QM occupation minus initial occupation) / "
            "RMS(estimate minus QM), per 100-au window and orbital"
        ),
        "strict_pass_rule": (
            "minimum Q over every active orbital and every window > 10"
        ),
        "cases": {},
    }
    rows: list[dict[str, object]] = []
    for label, estimate in estimates.items():
        result = metrics(qm, estimate)
        result["resources"] = resources(CASES[label].parent)
        result["strict_pass"] = result["Q_min_active_orbital_window"] > 10.0
        report["cases"][label] = result
        rows.append(
            {
                "case": label,
                "states": result["resources"]["reference_states"],
                "Q_global": result["Q_active_global"],
                "Q_min_active_window": result["Q_min_active_window"],
                "Q_min_active_orbital_window": result[
                    "Q_min_active_orbital_window"
                ],
                "max_abs_orbital_error": result["max_abs_orbital_error"],
                "reference_seconds": result["resources"]["reference_seconds"],
                "max_rank_peak_rss_kib": result["resources"][
                    "max_rank_peak_rss_kib"
                ],
                "strict_pass": result["strict_pass"],
            }
        )

    (ROOT / "validation.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    with (ROOT / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    time = qm[:, 0]
    initial = qm[0, 4 : 4 + NORB]
    colors = ["#d17a22", "#25867d", "#6e58a5"]
    fig, axes = plt.subplots(5, 2, figsize=(13, 14), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        axis.plot(
            time,
            qm[:, 4 + orbital] - initial[orbital],
            color="black",
            linewidth=1.15,
            label="Grid QM",
        )
        for color, (label, estimate) in zip(colors, estimates.items()):
            axis.plot(
                time,
                estimate[:, 4 + orbital] - initial[orbital],
                color=color,
                linewidth=0.7,
                alpha=0.85,
                label=label,
            )
        axis.set_title(f"Orbital {orbital}")
        axis.grid(alpha=0.18)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for axis in axes[-1]:
        axis.set_xlabel("Time (a.u.)")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(
        "10 orbitals / 5 electrons: bounded references vs grid QM",
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(ROOT / "reference_qm_all_orbitals.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for color, (label, estimate) in zip(colors, estimates.items()):
        result = report["cases"][label]
        centers = [
            0.5 * (window["start"] + window["stop"])
            for window in result["windows"]
        ]
        axes[0].semilogy(
            centers,
            [window["Q_weakest_active"] for window in result["windows"]],
            color=color,
            linewidth=1.1,
            label=label,
        )
        error = np.max(np.abs(estimate[:, 4:] - qm[:, 4:]), axis=1)
        axes[1].semilogy(
            time,
            np.maximum(error, 1.0e-18),
            color=color,
            linewidth=0.8,
            label=label,
        )
    axes[0].axhline(10.0, color="#9b3030", linestyle="--", linewidth=1.0)
    axes[0].set(
        xlabel="100-a.u. window midpoint",
        ylabel="Weakest active-orbital Q",
    )
    axes[1].set(
        xlabel="Time (a.u.)",
        ylabel="Largest absolute orbital error",
    )
    for axis in axes:
        axis.grid(alpha=0.18)
        axis.legend(fontsize=8)
    fig.suptitle("Strict 2500-a.u. validation; no smoothing or run averaging")
    fig.tight_layout()
    fig.savefig(ROOT / "reference_qm_quality.png", dpi=180)
    plt.close(fig)

    print(json.dumps(rows, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
