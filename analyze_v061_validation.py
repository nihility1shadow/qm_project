from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent / "variance_v061_validation_20260809"
GROUPS = {
    "v0.60 baseline": ["636252", "636253", "636254"],
    "v0.61 B=1": ["636255-v061-b1", "636256-v061-b1", "636257-v061-b1"],
    "v0.61 B=4": ["636258-v061-b4", "636259-v061-b4", "636260-v061-b4"],
    "v0.61 FB B=4": [
        "636263-v061fb-b4",
        "636264-v061fb-b4",
        "636265-v061fb-b4",
    ],
}
WINDOWS = ((0.0, 10.0), (10.0, 20.0), (20.0, 25.0), (0.0, 25.0))


def load_group(run_dirs: list[str]) -> tuple[np.ndarray, np.ndarray]:
    arrays = []
    for run_dir in run_dirs:
        path = ROOT / run_dir / "ahm-sepmb-s10-n5-100000.dat"
        arrays.append(np.loadtxt(path, comments="#"))
    times = arrays[0][:, 0]
    if any(not np.array_equal(times, array[:, 0]) for array in arrays[1:]):
        raise ValueError("time grids differ")
    return times, np.stack([array[:, 4:14] for array in arrays])


def summarize(name: str, times: np.ndarray, values: np.ndarray) -> list[dict[str, float | str]]:
    mean = values.mean(axis=0)
    repeat_std = values.std(axis=0, ddof=1)
    initial = np.r_[np.ones(5), np.zeros(5)]
    rows: list[dict[str, float | str]] = []
    for start, end in WINDOWS:
        mask = (times >= start) & (times <= end)
        noise_rms = float(np.sqrt(np.mean(repeat_std[mask] ** 2)))
        signal_rms = float(np.sqrt(np.mean((mean[mask] - initial) ** 2)))
        rows.append(
            {
                "group": name,
                "window": f"{start:g}-{end:g}",
                "noise_rms": noise_rms,
                "noise_mean_abs": float(np.mean(np.abs(repeat_std[mask]))),
                "noise_max": float(np.max(np.abs(repeat_std[mask]))),
                "signal_rms": signal_rms,
                "repeat_snr": signal_rms / noise_rms if noise_rms > 0 else float("inf"),
                "max_electron_sum_error": float(
                    np.max(np.abs(values[:, mask].sum(axis=2) - 5.0))
                ),
            }
        )
    return rows


def main() -> None:
    loaded = {name: load_group(run_dirs) for name, run_dirs in GROUPS.items()}
    qm_array = np.loadtxt(
        ROOT / "636261-qm" / "ahm-qm-s10-n5.dat", comments="#"
    )
    qm_times = qm_array[:, 0]
    qm_values = qm_array[:, 4:14]
    if any(not np.array_equal(times, qm_times) for times, _ in loaded.values()):
        raise ValueError("QM and Poisson time grids differ")
    rows = [
        row
        for name, (times, values) in loaded.items()
        for row in summarize(name, times, values)
    ]

    csv_path = ROOT / "validation_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / "validation_summary.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    colors = {
        "v0.60 baseline": "#555555",
        "v0.61 B=1": "#0072B2",
        "v0.61 B=4": "#D55E00",
        "v0.61 FB B=4": "#CC79A7",
    }
    fig, axes = plt.subplots(2, 5, figsize=(17, 7), sharex=True)
    initial = np.r_[np.ones(5), np.zeros(5)]
    for orbital, axis in enumerate(axes.flat):
        axis.plot(
            qm_times,
            qm_values[:, orbital] - initial[orbital],
            color="#009E73",
            lw=1.8,
            ls="--",
            label="exact QM",
        )
        for name, (times, values) in loaded.items():
            mean = values[:, :, orbital].mean(axis=0) - initial[orbital]
            std = values[:, :, orbital].std(axis=0, ddof=1)
            axis.plot(times, mean, color=colors[name], lw=1.35, label=name)
            axis.fill_between(times, mean - std, mean + std, color=colors[name], alpha=0.13)
        axis.axhline(0.0, color="#999999", lw=0.65)
        axis.set_title(f"orbital {orbital}")
        axis.grid(color="#dddddd", lw=0.55)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    for axis in axes[-1]:
        axis.set_xlabel("time (a.u.)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\Delta n_i(t)$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("s10n5 oxygen 2p: exact Bernoulli and stratified backward validation", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(ROOT / "all_orbitals_v060_v061.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(17, 7), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        axis.plot(
            qm_times,
            qm_values[:, orbital] - initial[orbital],
            color="#009E73",
            lw=1.8,
            ls="--",
            label="exact QM",
        )
        for name in ("v0.61 B=1", "v0.61 B=4", "v0.61 FB B=4"):
            times, values = loaded[name]
            mean = values[:, :, orbital].mean(axis=0) - initial[orbital]
            std = values[:, :, orbital].std(axis=0, ddof=1)
            axis.plot(times, mean, color=colors[name], lw=1.35, label=name)
            axis.fill_between(times, mean - std, mean + std, color=colors[name], alpha=0.13)
        axis.axhline(0.0, color="#999999", lw=0.65)
        axis.set_title(f"orbital {orbital}")
        axis.grid(color="#dddddd", lw=0.55)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    for axis in axes[-1]:
        axis.set_xlabel("time (a.u.)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\Delta n_i(t)$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("s10n5 oxygen 2p: v0.61 against exact QM", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(ROOT / "v061_vs_qm_all_orbitals.png", dpi=190)
    plt.close(fig)

    baseline_mean = loaded["v0.60 baseline"][1].mean(axis=0)
    comparisons = {}
    for name in GROUPS:
        mean = loaded[name][1].mean(axis=0)
        comparisons[name] = {
            "rms_mean_difference_from_v060": float(np.sqrt(np.mean((mean - baseline_mean) ** 2))),
            "max_mean_difference_from_v060": float(np.max(np.abs(mean - baseline_mean))),
            "rms_error_from_qm": float(np.sqrt(np.mean((mean - qm_values) ** 2))),
            "max_error_from_qm": float(np.max(np.abs(mean - qm_values))),
            "window_rms_error_from_qm": {
                f"{start:g}-{end:g}": float(
                    np.sqrt(np.mean((mean[(qm_times >= start) & (qm_times <= end)] -
                                     qm_values[(qm_times >= start) & (qm_times <= end)]) ** 2))
                )
                for start, end in WINDOWS
            },
        }
    (ROOT / "mean_change_from_v060.json").write_text(
        json.dumps(comparisons, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
