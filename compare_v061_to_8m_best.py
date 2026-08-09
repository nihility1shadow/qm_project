from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "compare_v061_vs_8m_best_20260809"
OLD_FILES = [
    ROOT / "scan_cloud_8m_t150_20260809" / "round4" / "cloud-runs" / job
    / "ahm-sepmb-s10-n5-8000000.dat"
    for job in ("636039", "636040", "636041")
] + [
    ROOT / "scan_cloud_8m_t150_20260809" / "round5" / "cloud-runs" / job
    / "ahm-sepmb-s10-n5-8000000.dat"
    for job in ("636076", "636077", "636078")
]
NEW_FILE = (
    OUT / "636268-v061-8m-best" / "ahm-sepmb-s10-n5-8000000.dat"
)
REPRESENTATIVE_OLD_JOB = "636040"
WINDOWS = tuple((float(start), float(start + 25)) for start in range(0, 150, 25))


def load(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.shape[1] != 14:
        raise ValueError(f"unexpected columns in {path}: {data.shape}")
    if not np.isfinite(data).all():
        raise ValueError(f"non-finite values in {path}")
    return data


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    old_data = [load(path) for path in OLD_FILES]
    new_data = load(NEW_FILE)
    times = new_data[:, 0]
    if any(not np.array_equal(times, data[:, 0]) for data in old_data):
        raise ValueError("old and new time grids differ")

    old_occ = np.stack([data[:, 4:14] for data in old_data])
    new_occ = new_data[:, 4:14]
    old_mean = old_occ.mean(axis=0)
    old_std = old_occ.std(axis=0, ddof=1)
    representative_index = next(
        index
        for index, path in enumerate(OLD_FILES)
        if path.parent.name == REPRESENTATIVE_OLD_JOB
    )
    old_representative = old_occ[representative_index]
    initial = np.r_[np.ones(5), np.zeros(5)]

    checks = {
        "parameters": {
            "wc_eV": 0.12,
            "eta": 2.8e-5,
            "delE_eV": -2.5,
            "Norb": 10,
            "Nel": 5,
            "ntraj_per_run": 8_000_000,
            "dt": 0.5,
            "tmax": 150.0,
        },
        "old_jobs": [636039, 636040, 636041, 636076, 636077, 636078],
        "new_job": 636268,
        "representative_old_job": int(REPRESENTATIVE_OLD_JOB),
        "representative_selection": "old 8M run with the lowest RMS distance to the six-run old mean",
        "old_repeats": 6,
        "new_repeats": 1,
        "new_max_particle_number_error": float(
            np.max(np.abs(new_occ.sum(axis=1) - 5.0))
        ),
        "old_max_particle_number_error": float(
            np.max(np.abs(old_occ.sum(axis=2) - 5.0))
        ),
        "new_all_finite": bool(np.isfinite(new_data).all()),
    }
    (OUT / "comparison_config_and_checks.json").write_text(
        json.dumps(checks, indent=2), encoding="utf-8"
    )

    rows = []
    for start, end in WINDOWS:
        mask = (times >= start) & (times <= end)
        difference = new_occ[mask] - old_mean[mask]
        same_8m_difference = new_occ[mask] - old_representative[mask]
        old_noise = old_std[mask]
        rows.append(
            {
                "window": f"{start:g}-{end:g}",
                "rms_new_8m_minus_old_8m": float(
                    np.sqrt(np.mean(same_8m_difference**2))
                ),
                "max_abs_new_8m_minus_old_8m": float(
                    np.max(np.abs(same_8m_difference))
                ),
                "rms_new_minus_old_mean": float(np.sqrt(np.mean(difference**2))),
                "max_abs_new_minus_old_mean": float(np.max(np.abs(difference))),
                "old_repeat_noise_rms": float(np.sqrt(np.mean(old_noise**2))),
                "rms_difference_over_old_noise": float(
                    np.sqrt(np.mean(difference**2))
                    / np.sqrt(np.mean(old_noise**2))
                    if np.any(old_noise) else float("inf")
                ),
            }
        )
    with (OUT / "window_comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "window_comparison.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8"
    )

    per_old = []
    for path, values in zip(OLD_FILES, old_occ):
        per_old.append(
            {
                "old_job": path.parent.name,
                "rms_new_minus_old_single": float(
                    np.sqrt(np.mean((new_occ - values) ** 2))
                ),
                "max_abs_new_minus_old_single": float(
                    np.max(np.abs(new_occ - values))
                ),
            }
        )
    with (OUT / "new_vs_each_old_run.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_old[0]))
        writer.writeheader()
        writer.writerows(per_old)

    def plot_same_sample_comparison(limit: float | None, filename: str) -> None:
        mask = times <= limit if limit is not None else np.ones(times.shape, dtype=bool)
        plot_times = times[mask]
        fig, axes = plt.subplots(2, 5, figsize=(18, 7.5), sharex=True)
        for orbital, axis in enumerate(axes.flat):
            old_single_delta = old_representative[mask, orbital] - initial[orbital]
            new_delta = new_occ[mask, orbital] - initial[orbital]
            old_mean_delta = old_mean[mask, orbital] - initial[orbital]
            axis.plot(
                plot_times,
                old_single_delta,
                color="#0072B2",
                lw=1.25,
                label="v0.60 old, 8M (job 636040)" if orbital == 0 else None,
            )
            axis.plot(
                plot_times,
                new_delta,
                color="#D55E00",
                lw=1.35,
                label="v0.61 fixed, 8M (job 636268)" if orbital == 0 else None,
            )
            axis.plot(
                plot_times,
                old_mean_delta,
                color="#555555",
                lw=0.95,
                ls="--",
                alpha=0.72,
                label="v0.60 mean of 6 runs (reference)" if orbital == 0 else None,
            )
            axis.axhline(0.0, color="#bdbdbd", lw=0.6)
            axis.set_title(f"orbital {orbital}")
            axis.grid(color="#e0e0e0", lw=0.5)
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        for axis in axes[-1]:
            axis.set_xlabel("time (a.u.)")
        for axis in axes[:, 0]:
            axis.set_ylabel(r"$\Delta n_i(t)$")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.945),
            ncol=3,
            frameon=False,
        )
        interval = f"0-{limit:g} a.u." if limit is not None else "0-150 a.u."
        fig.suptitle(
            f"Same-sample comparison ({interval}): 8M trajectories per solid curve\n"
            "wc=0.12 eV, eta=2.8e-5, delE=-2.5 eV, Norb=10, Nel=5",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.855))
        fig.savefig(OUT / filename, dpi=190)
        plt.close(fig)

    plot_same_sample_comparison(
        None, "same_8m_old636040_vs_new636268_full_t150.png"
    )
    plot_same_sample_comparison(
        60.0, "same_8m_old636040_vs_new636268_zoom_t60.png"
    )

    fig, axes = plt.subplots(2, 5, figsize=(18, 7.5), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        for repeat in old_occ:
            axis.plot(
                times,
                repeat[:, orbital] - initial[orbital],
                color="#9a9a9a",
                lw=0.65,
                alpha=0.34,
            )
        old_delta = old_mean[:, orbital] - initial[orbital]
        new_delta = new_occ[:, orbital] - initial[orbital]
        axis.fill_between(
            times,
            old_delta - old_std[:, orbital],
            old_delta + old_std[:, orbital],
            color="#777777",
            alpha=0.14,
            label="v0.60 old +/- 1 SD" if orbital == 0 else None,
        )
        axis.plot(
            times,
            old_delta,
            color="#333333",
            lw=1.25,
            ls="--",
            label="v0.60 mean (6 runs)" if orbital == 0 else None,
        )
        axis.plot(
            times,
            new_delta,
            color="#D55E00",
            lw=1.45,
            label="v0.61 single run" if orbital == 0 else None,
        )
        axis.axhline(0.0, color="#bdbdbd", lw=0.6)
        axis.set_title(f"orbital {orbital}")
        axis.grid(color="#e0e0e0", lw=0.5)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    for axis in axes[-1]:
        axis.set_xlabel("time (a.u.)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\Delta n_i(t)$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        "8M trajectories, wc=0.12 eV, eta=2.8e-5, delE=-2.5 eV",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(OUT / "all_orbitals_v060_best_vs_v061_8m.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(18, 7.5), sharex=True)
    difference = new_occ - old_mean
    for orbital, axis in enumerate(axes.flat):
        axis.plot(times, difference[:, orbital], color="#0072B2", lw=1.3)
        axis.fill_between(
            times,
            -old_std[:, orbital],
            old_std[:, orbital],
            color="#777777",
            alpha=0.18,
        )
        axis.axhline(0.0, color="#555555", lw=0.65)
        axis.set_title(f"orbital {orbital}")
        axis.grid(color="#e0e0e0", lw=0.5)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    for axis in axes[-1]:
        axis.set_xlabel("time (a.u.)")
    for axis in axes[:, 0]:
        axis.set_ylabel("v0.61 - v0.60 mean")
    fig.suptitle("Difference; gray band is old +/- 1 repeat SD", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / "all_orbitals_difference_vs_old_noise.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    main()
