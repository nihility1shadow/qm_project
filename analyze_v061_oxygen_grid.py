from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SCAN = ROOT / "scan_v061_oxygen_grid_20260809"
INITIAL = np.r_[np.ones(5), np.zeros(5)]


def load_data(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] != 14:
        raise ValueError(f"unexpected data shape in {path}: {data.shape}")
    if not np.isfinite(data).all():
        raise ValueError(f"non-finite value in {path}")
    return data


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def data_path(row: dict[str, str], stage: str) -> Path:
    repeat = "exact" if row["method"] == "qm" else f"rep{row['repeat']}"
    filename = (
        "ahm-qm-s10-n5.dat"
        if row["method"] == "qm"
        else f"ahm-sepmb-s10-n5-{row['ntraj']}.dat"
    )
    return SCAN / f"{stage}_data" / row["case_id"] / row["method"] / repeat / filename


def analyze(stage: str = "coarse") -> None:
    with (SCAN / f"{stage}_submission.csv").open(encoding="utf-8-sig") as handle:
        submissions = list(csv.DictReader(handle))

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in submissions:
        grouped.setdefault(row["case_id"], []).append(row)

    metrics: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []
    cached: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

    for case_id, rows in grouped.items():
        poisson_rows = sorted(
            (row for row in rows if row["method"] == "v061"),
            key=lambda row: int(row["repeat"]),
        )
        qm_rows = [row for row in rows if row["method"] == "qm"]
        if len(poisson_rows) < 2 or len(qm_rows) != 1:
            continue

        poisson_data = [load_data(data_path(row, stage)) for row in poisson_rows]
        qm_data = load_data(data_path(qm_rows[0], stage))
        times = qm_data[:, 0]
        if any(not np.array_equal(times, values[:, 0]) for values in poisson_data):
            raise ValueError(f"time-grid mismatch for {case_id}")

        poisson_occ = np.stack([values[:, 4:14] for values in poisson_data])
        qm_occ = qm_data[:, 4:14]
        poisson_mean = poisson_occ.mean(axis=0)
        poisson_std = poisson_occ.std(axis=0, ddof=1)
        mean_noise = poisson_std / np.sqrt(len(poisson_data))
        signal = qm_occ - INITIAL
        error = poisson_mean - qm_occ

        signal_rms = rms(signal)
        signal_peak = float(np.max(np.abs(signal)))
        error_rms = rms(error)
        noise_rms = rms(mean_noise)
        q_accuracy = signal_rms / max(error_rms, 1e-30)
        q_repeat = signal_rms / max(noise_rms, 1e-30)
        conservative_q = min(q_accuracy, q_repeat)

        base = poisson_rows[0]
        metrics.append(
            {
                "stage": stage,
                "case_id": case_id,
                "version": "v0.61",
                "wc_eV": float(base["wc_eV"]),
                "eta": float(base["eta"]),
                "delE_eV": float(base["delE_eV"]),
                "ntraj_per_repeat": int(base["ntraj"]),
                "independent_repeats": len(poisson_data),
                "tmax": float(times[-1]),
                "qm_signal_rms": signal_rms,
                "qm_signal_peak": signal_peak,
                "poisson_rms_vs_qm": error_rms,
                "repeat_noise_of_mean_rms": noise_rms,
                "Q_accuracy": q_accuracy,
                "Q_repeat": q_repeat,
                "Q_conservative": conservative_q,
                "projected_ntraj_Q10": int(
                    np.ceil(int(base["ntraj"]) * (10.0 / conservative_q) ** 2)
                ),
                "occupation_min": float(poisson_occ.min()),
                "occupation_max": float(poisson_occ.max()),
                "max_particle_error": float(
                    np.max(np.abs(poisson_occ.sum(axis=2) - 5.0))
                ),
            }
        )

        windows = tuple(
            (float(start), min(float(start + 10), float(times[-1])))
            for start in range(0, int(np.ceil(times[-1])), 10)
        )
        for start, end in windows:
            mask = (times >= start) & (times <= end)
            window_signal = rms(signal[mask])
            window_error = rms(error[mask])
            window_noise = rms(mean_noise[mask])
            window_rows.append(
                {
                    "stage": stage,
                    "case_id": case_id,
                    "window": f"{start:g}-{end:g}",
                    "wc_eV": float(base["wc_eV"]),
                    "eta": float(base["eta"]),
                    "signal_rms": window_signal,
                    "error_rms_vs_qm": window_error,
                    "repeat_noise_rms": window_noise,
                    "Q_accuracy": window_signal / max(window_error, 1e-30),
                    "Q_repeat": window_signal / max(window_noise, 1e-30),
                    "Q_conservative": min(
                        window_signal / max(window_error, 1e-30),
                        window_signal / max(window_noise, 1e-30),
                    ),
                }
            )

        cached[case_id] = (times, poisson_mean, poisson_std, qm_occ)

    if not metrics:
        raise RuntimeError(f"no complete {stage} cases found")

    max_signal = max(float(row["qm_signal_rms"]) for row in metrics)
    windows_by_case: dict[str, list[dict[str, object]]] = {}
    for row in window_rows:
        windows_by_case.setdefault(str(row["case_id"]), []).append(row)
    for row in metrics:
        case_windows = windows_by_case[str(row["case_id"])]
        min_window_q = min(float(item["Q_conservative"]) for item in case_windows)
        row["min_window_Q"] = min_window_q
        row["selection_score"] = float(row["Q_conservative"]) * np.sqrt(
            float(row["qm_signal_rms"]) / max(max_signal, 1e-30)
        )

    metrics.sort(key=lambda row: float(row["selection_score"]), reverse=True)
    for rank, row in enumerate(metrics, start=1):
        row["rank"] = rank

    field_order = ["rank"] + [key for key in metrics[0] if key != "rank"]
    metrics = [{key: row[key] for key in field_order} for row in metrics]
    write_rows(SCAN / f"{stage}_combination_metrics.csv", metrics)
    write_rows(SCAN / f"{stage}_window_metrics.csv", window_rows)
    (SCAN / f"{stage}_combination_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )

    wc_values = sorted({float(row["wc_eV"]) for row in metrics})
    eta_values = sorted({float(row["eta"]) for row in metrics})
    score = np.full((len(eta_values), len(wc_values)), np.nan)
    quality = np.full_like(score, np.nan)
    signal_map = np.full_like(score, np.nan)
    for row in metrics:
        i = eta_values.index(float(row["eta"]))
        j = wc_values.index(float(row["wc_eV"]))
        score[i, j] = float(row["selection_score"])
        quality[i, j] = float(row["Q_conservative"])
        signal_map[i, j] = float(row["qm_signal_rms"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.3))
    for axis, values, title, colorbar_label in (
        (axes[0], score, "selection score", "Q weighted by signal"),
        (axes[1], quality, "conservative Q", "min(accuracy Q, repeat Q)"),
        (axes[2], signal_map, "exact-QM signal", "RMS occupation change"),
    ):
        image = axis.imshow(values, origin="lower", aspect="auto", cmap="viridis")
        axis.set_xticks(range(len(wc_values)), [f"{value:g}" for value in wc_values])
        axis.set_yticks(
            range(len(eta_values)), [f"{value:.3g}" for value in eta_values]
        )
        axis.set_xlabel("wc (eV)")
        axis.set_ylabel("eta")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, label=colorbar_label)
    first = metrics[0]
    fig.suptitle(
        f"v0.61 oxygen 2p {stage} grid: delE=-17.1948749688 eV, "
        f"{first['independent_repeats']} x {int(first['ntraj_per_repeat']):,} "
        f"trajectories, t=0-{first['tmax']:g} a.u."
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(SCAN / f"{stage}_parameter_maps.png", dpi=190)
    plt.close(fig)

    plots = SCAN / f"{stage}_top_plots"
    plots.mkdir(exist_ok=True)
    for row in metrics[:8]:
        case_id = str(row["case_id"])
        times, mean, std, qm = cached[case_id]
        fig, axes = plt.subplots(2, 5, figsize=(18, 7.3), sharex=True)
        for orbital, axis in enumerate(axes.flat):
            offset = INITIAL[orbital]
            axis.fill_between(
                times,
                mean[:, orbital] - offset - std[:, orbital],
                mean[:, orbital] - offset + std[:, orbital],
                color="#D55E00",
                alpha=0.15,
            )
            axis.plot(
                times,
                mean[:, orbital] - offset,
                color="#D55E00",
                lw=1.2,
                label="v0.61 mean" if orbital == 0 else None,
            )
            axis.plot(
                times,
                qm[:, orbital] - offset,
                color="#009E73",
                lw=1.45,
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
            ncol=2,
            frameon=False,
        )
        fig.suptitle(
            f"rank {row['rank']}: wc={row['wc_eV']:g} eV, "
            f"eta={row['eta']:.3g}, Q={row['Q_conservative']:.3g}",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.865))
        fig.savefig(plots / f"rank{int(row['rank']):02d}_{case_id}.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    analyze(sys.argv[1] if len(sys.argv) > 1 else "coarse")
