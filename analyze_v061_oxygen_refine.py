from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SCAN = ROOT / "scan_v061_oxygen_local_refine_20260810"
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


def quality(signal: np.ndarray, error: np.ndarray, noise: np.ndarray) -> tuple[float, float, float]:
    signal_rms = rms(signal)
    q_accuracy = signal_rms / max(rms(error), 1e-30)
    q_repeat = signal_rms / max(rms(noise), 1e-30)
    return q_accuracy, q_repeat, min(q_accuracy, q_repeat)


def analyze(stage: str) -> None:
    submission_path = SCAN / f"{stage}_submission.csv"
    with submission_path.open(encoding="utf-8-sig") as handle:
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
        q_accuracy, q_repeat, conservative_q = quality(signal, error, mean_noise)

        interval_quality: dict[str, float] = {}
        for label, start, end in (
            ("Q_0_40", 0.0, 40.0),
            ("Q_40_end", 40.0, float(times[-1])),
            ("Q_0_60", 0.0, min(60.0, float(times[-1]))),
        ):
            mask = (times >= start) & (times <= end)
            interval_quality[label] = quality(signal[mask], error[mask], mean_noise[mask])[2]

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
                "qm_signal_rms": rms(signal),
                "qm_signal_peak": float(np.max(np.abs(signal))),
                "poisson_rms_vs_qm": rms(error),
                "repeat_noise_of_mean_rms": rms(mean_noise),
                "Q_accuracy": q_accuracy,
                "Q_repeat": q_repeat,
                "Q_conservative": conservative_q,
                **interval_quality,
                "projected_ntraj_Q10": int(
                    math.ceil(int(base["ntraj"]) * (10.0 / conservative_q) ** 2)
                ),
                "occupation_min": float(poisson_occ.min()),
                "occupation_max": float(poisson_occ.max()),
                "max_particle_error": float(
                    np.max(np.abs(poisson_occ.sum(axis=2) - 5.0))
                ),
            }
        )

        for start in range(0, int(math.ceil(times[-1])), 10):
            end = min(float(start + 10), float(times[-1]))
            mask = (times >= start) & (times <= end)
            qa, qr, qc = quality(signal[mask], error[mask], mean_noise[mask])
            window_rows.append(
                {
                    "stage": stage,
                    "case_id": case_id,
                    "window": f"{start:g}-{end:g}",
                    "wc_eV": float(base["wc_eV"]),
                    "eta": float(base["eta"]),
                    "signal_rms": rms(signal[mask]),
                    "error_rms_vs_qm": rms(error[mask]),
                    "repeat_noise_rms": rms(mean_noise[mask]),
                    "Q_accuracy": qa,
                    "Q_repeat": qr,
                    "Q_conservative": qc,
                }
            )

        cached[case_id] = (times, poisson_mean, poisson_std, qm_occ)

    if not metrics:
        raise RuntimeError(f"no complete {stage} cases found")

    windows_by_case: dict[str, list[dict[str, object]]] = {}
    for row in window_rows:
        windows_by_case.setdefault(str(row["case_id"]), []).append(row)
    max_signal = max(float(row["qm_signal_rms"]) for row in metrics)
    for row in metrics:
        row["min_window_Q"] = min(
            float(item["Q_conservative"])
            for item in windows_by_case[str(row["case_id"])]
        )
        robust_q = min(float(row["Q_conservative"]), float(row["Q_0_40"]))
        row["selection_score"] = robust_q * math.sqrt(
            float(row["qm_signal_rms"]) / max(max_signal, 1e-30)
        )

    metrics.sort(key=lambda row: float(row["selection_score"]), reverse=True)
    for rank, row in enumerate(metrics, start=1):
        row["rank"] = rank
    field_order = ["rank"] + [key for key in metrics[0] if key != "rank"]
    metrics = [{key: row[key] for key in field_order} for row in metrics]

    retained_count = max(1, math.ceil(len(metrics) * 0.30))
    retained = [dict(row, selection_status="retained_top30") for row in metrics[:retained_count]]
    write_rows(SCAN / f"{stage}_metrics.csv", metrics)
    write_rows(SCAN / f"{stage}_window_metrics.csv", window_rows)
    write_rows(SCAN / f"{stage}_retained_top30.csv", retained)
    (SCAN / f"{stage}_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )

    wc_values = sorted({float(row["wc_eV"]) for row in metrics})
    eta_values = sorted({float(row["eta"]) for row in metrics})
    maps = {
        "selection score": np.full((len(eta_values), len(wc_values)), np.nan),
        "conservative Q": np.full((len(eta_values), len(wc_values)), np.nan),
        "Q, 40-end": np.full((len(eta_values), len(wc_values)), np.nan),
        "projected paths for Q=10": np.full((len(eta_values), len(wc_values)), np.nan),
    }
    for row in metrics:
        i = eta_values.index(float(row["eta"]))
        j = wc_values.index(float(row["wc_eV"]))
        maps["selection score"][i, j] = float(row["selection_score"])
        maps["conservative Q"][i, j] = float(row["Q_conservative"])
        maps["Q, 40-end"][i, j] = float(row["Q_40_end"])
        maps["projected paths for Q=10"][i, j] = float(row["projected_ntraj_Q10"])

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for axis, (title, values) in zip(axes.flat, maps.items()):
        image = axis.imshow(values, origin="lower", aspect="auto", cmap="viridis")
        axis.set_xticks(range(len(wc_values)), [f"{value:g}" for value in wc_values])
        axis.set_yticks(range(len(eta_values)), [f"{value:.2g}" for value in eta_values])
        axis.set_xlabel("wc (eV)")
        axis.set_ylabel("eta")
        axis.set_title(title)
        fig.colorbar(image, ax=axis)
    first = metrics[0]
    fig.suptitle(
        f"v0.61 local {stage}: delE=-17.1948749688 eV, "
        f"{first['independent_repeats']} x {int(first['ntraj_per_repeat']):,} paths, "
        f"t=0-{first['tmax']:g} a.u."
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(SCAN / f"{stage}_parameter_maps.png", dpi=190)
    plt.close(fig)

    plot_root = SCAN / f"{stage}_top_plots"
    plot_root.mkdir(exist_ok=True)
    for row in metrics[:retained_count]:
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
                alpha=0.14,
            )
            axis.plot(times, mean[:, orbital] - offset, color="#D55E00", lw=1.1)
            axis.plot(times, qm[:, orbital] - offset, color="#0072B2", lw=1.35)
            axis.axhline(0.0, color="#b5b5b5", lw=0.6)
            axis.grid(color="#e2e2e2", lw=0.5)
            axis.set_title(f"orbital {orbital}")
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        for axis in axes[-1]:
            axis.set_xlabel("time (a.u.)")
        for axis in axes[:, 0]:
            axis.set_ylabel(r"$\Delta n_i(t)$")
        fig.suptitle(
            f"rank {row['rank']}: wc={row['wc_eV']:g} eV, eta={row['eta']:.3g}, "
            f"Q={row['Q_conservative']:.3g}, Q(40-end)={row['Q_40_end']:.3g}",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(plot_root / f"rank{int(row['rank']):02d}_{case_id}.png", dpi=175)
        plt.close(fig)

    print(
        json.dumps(
            {
                "stage": stage,
                "combinations": len(metrics),
                "retained": retained_count,
                "best": metrics[0],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    analyze(sys.argv[1] if len(sys.argv) > 1 else "screen")
