from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SCAN = ROOT / "scan_v061_oxygen_grid_20260809"
INITIAL = np.r_[np.ones(5), np.zeros(5)]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load(path: Path) -> np.ndarray:
    values = np.loadtxt(path, comments="#")
    if values.shape[1] != 14 or not np.isfinite(values).all():
        raise ValueError(f"invalid data in {path}")
    return values


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def reference_metrics(
    version: str,
    arrays: list[np.ndarray],
    qm: np.ndarray,
    ntraj: int,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    times = arrays[0][:, 0]
    if any(not np.array_equal(times, values[:, 0]) for values in arrays):
        raise ValueError(f"reference time-grid mismatch for {version}")
    qm_window = qm[np.isin(qm[:, 0], times)]
    if not np.array_equal(times, qm_window[:, 0]):
        raise ValueError(f"QM time-grid mismatch for {version}")
    occupations = np.stack([values[:, 4:14] for values in arrays])
    mean = occupations.mean(axis=0)
    std = occupations.std(axis=0, ddof=1)
    noise_of_mean = std / np.sqrt(len(arrays))
    qm_occ = qm_window[:, 4:14]
    signal = qm_occ - INITIAL
    error = mean - qm_occ
    signal_rms = rms(signal)
    error_rms = rms(error)
    noise_rms = rms(noise_of_mean)
    q_accuracy = signal_rms / max(error_rms, 1e-30)
    q_repeat = signal_rms / max(noise_rms, 1e-30)
    q = min(q_accuracy, q_repeat)
    return (
        {
            "catalog_id": f"reference25_{version}",
            "stage": "old_parameter_reference_0_25",
            "selection_status": "reference",
            "version": version,
            "case_id": "wc4_eta0p00013758",
            "wc_eV": 4.0,
            "eta": 0.00013758,
            "delE_eV": -17.194874968839816,
            "ntraj_per_repeat": ntraj,
            "independent_repeats": len(arrays),
            "tmax": float(times[-1]),
            "qm_signal_rms": signal_rms,
            "qm_signal_peak": float(np.max(np.abs(signal))),
            "poisson_rms_vs_qm": error_rms,
            "repeat_noise_of_mean_rms": noise_rms,
            "Q_accuracy": q_accuracy,
            "Q_repeat": q_repeat,
            "Q_conservative": q,
            "projected_ntraj_Q10": int(math.ceil(ntraj * (10.0 / q) ** 2)),
            "occupation_min": float(occupations.min()),
            "occupation_max": float(occupations.max()),
            "max_particle_error": float(
                np.max(np.abs(occupations.sum(axis=2) - 5.0))
            ),
            "min_window_Q": "",
            "selection_score": "",
            "rank": "",
            "source": "direct old-parameter comparison against matching QM",
        },
        mean,
        std,
    )


def normalized_grid_row(
    row: dict[str, str], retained_count: int
) -> dict[str, object]:
    rank = int(row["rank"])
    return {
        "catalog_id": f"{row['stage']}_{row['case_id']}",
        "stage": row["stage"],
        "selection_status": "retained_top30" if rank <= retained_count else "screened_out",
        "version": row["version"],
        "case_id": row["case_id"],
        "wc_eV": float(row["wc_eV"]),
        "eta": float(row["eta"]),
        "delE_eV": float(row["delE_eV"]),
        "ntraj_per_repeat": int(row["ntraj_per_repeat"]),
        "independent_repeats": int(row["independent_repeats"]),
        "tmax": float(row["tmax"]),
        "qm_signal_rms": float(row["qm_signal_rms"]),
        "qm_signal_peak": float(row["qm_signal_peak"]),
        "poisson_rms_vs_qm": float(row["poisson_rms_vs_qm"]),
        "repeat_noise_of_mean_rms": float(row["repeat_noise_of_mean_rms"]),
        "Q_accuracy": float(row["Q_accuracy"]),
        "Q_repeat": float(row["Q_repeat"]),
        "Q_conservative": float(row["Q_conservative"]),
        "projected_ntraj_Q10": int(row["projected_ntraj_Q10"]),
        "occupation_min": float(row["occupation_min"]),
        "occupation_max": float(row["occupation_max"]),
        "max_particle_error": float(row["max_particle_error"]),
        "min_window_Q": float(row["min_window_Q"]),
        "selection_score": float(row["selection_score"]),
        "rank": rank,
        "source": f"v0.61 {row['stage']} grid with matching exact QM",
    }


def main() -> None:
    coarse = read_csv(SCAN / "coarse_combination_metrics.csv")
    fine = read_csv(SCAN / "fine_combination_metrics.csv")
    coarse_keep = math.ceil(0.30 * len(coarse))
    fine_keep = math.ceil(0.30 * len(fine))

    catalog = [normalized_grid_row(row, coarse_keep) for row in coarse]
    catalog.extend(normalized_grid_row(row, fine_keep) for row in fine)

    old_paths = [
        ROOT
        / "scan_o2p_lit_t25_n1000000_20260713"
        / f"D1_w4_e1376_rep{repeat}"
        / "ahm-sepmb-s10-n5-1000000.dat"
        for repeat in (1, 2)
    ]
    reference_root = (
        SCAN
        / "coarse_data"
        / "ref_old_v060_wc4_eta0p00013758"
    )
    new_paths = [
        reference_root / "v061" / f"rep{repeat}" / "ahm-sepmb-s10-n5-500000.dat"
        for repeat in (1, 2)
    ]
    qm = load(reference_root / "qm" / "exact" / "ahm-qm-s10-n5.dat")
    old_arrays = [load(path) for path in old_paths]
    new_arrays = [load(path) for path in new_paths]
    new_arrays = [values[values[:, 0] <= 25.0] for values in new_arrays]

    old_row, old_mean, old_std = reference_metrics(
        "v0.60", old_arrays, qm, 1_000_000
    )
    new_row, new_mean, new_std = reference_metrics(
        "v0.61", new_arrays, qm, 500_000
    )
    catalog.extend([old_row, new_row])

    write_csv(SCAN / "oxygen_grid_combination_catalog.csv", catalog)
    write_csv(SCAN / "old_parameter_direct_comparison.csv", [old_row, new_row])

    retained_fine = [
        row for row in catalog
        if row["stage"] == "fine" and row["selection_status"] == "retained_top30"
    ]
    write_csv(SCAN / "fine_retained_top30.csv", retained_fine)

    qm_25 = qm[qm[:, 0] <= 25.0]
    times = qm_25[:, 0]
    fig, axes = plt.subplots(2, 5, figsize=(18, 7.5), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        offset = INITIAL[orbital]
        axis.fill_between(
            times,
            old_mean[:, orbital] - offset - old_std[:, orbital],
            old_mean[:, orbital] - offset + old_std[:, orbital],
            color="#0072B2",
            alpha=0.12,
        )
        axis.fill_between(
            times,
            new_mean[:, orbital] - offset - new_std[:, orbital],
            new_mean[:, orbital] - offset + new_std[:, orbital],
            color="#D55E00",
            alpha=0.12,
        )
        axis.plot(
            times,
            old_mean[:, orbital] - offset,
            color="#0072B2",
            lw=1.15,
            label="v0.60: 2 x 1M" if orbital == 0 else None,
        )
        axis.plot(
            times,
            new_mean[:, orbital] - offset,
            color="#D55E00",
            lw=1.25,
            label="v0.61: 2 x 500k" if orbital == 0 else None,
        )
        axis.plot(
            times,
            qm_25[:, 4 + orbital] - offset,
            color="#009E73",
            lw=1.5,
            label="exact QM" if orbital == 0 else None,
        )
        axis.axhline(0.0, color="#bbbbbb", lw=0.6)
        axis.grid(color="#e1e1e1", lw=0.5)
        axis.set_title(f"orbital {orbital}")
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
        "Old-reference parameters: wc=4 eV, eta=1.3758e-4, "
        "delE=-17.1948749688 eV, t=0-25 a.u.",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    fig.savefig(SCAN / "old_v060_v061_qm_reference_t25.png", dpi=190)
    plt.close(fig)

    labels = [f"{row['wc_eV']:g}, {float(row['eta']):.1e}" for row in retained_fine]
    q_values = [float(row["Q_conservative"]) for row in retained_fine]
    projected = [int(row["projected_ntraj_Q10"]) / 1e6 for row in retained_fine]
    x = np.arange(len(labels))
    fig, axis_q = plt.subplots(figsize=(13, 6.2))
    axis_q.bar(x, q_values, color="#0072B2", alpha=0.82)
    axis_q.set_ylabel("conservative Q at 3 x 1M, t=0-75")
    axis_q.set_xticks(x, labels, rotation=28, ha="right")
    axis_q.set_xlabel("wc (eV), eta")
    axis_q.grid(axis="y", color="#dddddd", lw=0.6)
    axis_n = axis_q.twinx()
    axis_n.plot(x, projected, color="#D55E00", marker="o", lw=1.5)
    axis_n.set_ylabel("projected trajectories per repeat for Q=10 (million)")
    axis_q.set_title("Fine-grid retained top 30%")
    fig.tight_layout()
    fig.savefig(SCAN / "fine_retained_top30_summary.png", dpi=190)
    plt.close(fig)

    summary = {
        "fixed_delE_eV": -17.194874968839816,
        "coarse_combinations": len(coarse),
        "fine_combinations": len(fine),
        "catalog_rows": len(catalog),
        "coarse_retained": coarse_keep,
        "fine_retained": fine_keep,
        "best_fine": retained_fine[0],
        "old_reference": old_row,
        "new_reference": new_row,
    }
    (SCAN / "oxygen_grid_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
