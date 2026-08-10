from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SCAN = ROOT / "scan_v061_oxygen_t200_q10_20260810"
INITIAL = np.r_[np.ones(5), np.zeros(5)]


def load_data(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] != 14:
        raise ValueError(f"unexpected data shape in {path}: {data.shape}")
    if not np.isfinite(data).all():
        raise ValueError(f"non-finite data in {path}")
    return data


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(signal: np.ndarray, error: np.ndarray, noise: np.ndarray) -> tuple[float, float, float]:
    signal_rms = rms(signal)
    q_accuracy = signal_rms / max(rms(error), 1e-30)
    q_repeat = signal_rms / max(rms(noise), 1e-30)
    return q_accuracy, q_repeat, min(q_accuracy, q_repeat)


def interval_quality(
    times: np.ndarray,
    signal: np.ndarray,
    error: np.ndarray,
    noise: np.ndarray,
    start: float,
    end: float,
) -> tuple[float, float, float]:
    mask = (times >= start) & (times <= end)
    if not np.any(mask):
        return math.nan, math.nan, math.nan
    return quality(signal[mask], error[mask], noise[mask])


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(stage: str) -> None:
    submission = SCAN / f"{stage}_submission.csv"
    with submission.open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    poisson_rows: dict[str, list[dict[str, str]]] = {}
    qm_rows: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["method"] == "qm":
            qm_rows[row["param_id"]] = row
        else:
            poisson_rows.setdefault(row["case_id"], []).append(row)

    metrics: list[dict[str, object]] = []
    window_rows: list[dict[str, object]] = []
    cached: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for case_id, case_rows in poisson_rows.items():
        case_rows.sort(key=lambda row: int(row["repeat"]))
        if len(case_rows) < 2:
            continue
        base = case_rows[0]
        param_id = base["param_id"]
        qm_path = SCAN / f"{stage}_data" / "qm" / param_id / "ahm-qm-s10-n5.dat"
        qm_data = load_data(qm_path)
        poisson_data = [
            load_data(
                SCAN
                / f"{stage}_data"
                / case_id
                / "v061"
                / f"rep{row['repeat']}"
                / f"ahm-sepmb-s10-n5-{row['ntraj']}.dat"
            )
            for row in case_rows
        ]

        times = qm_data[:, 0]
        if any(not np.array_equal(times, data[:, 0]) for data in poisson_data):
            raise ValueError(f"time grid mismatch: {case_id}")
        if not math.isclose(float(times[-1]), 200.0, abs_tol=1e-9):
            raise ValueError(f"expected tmax=200 for {case_id}, found {times[-1]}")

        occupations = np.stack([data[:, 4:14] for data in poisson_data])
        qm = qm_data[:, 4:14]
        mean = occupations.mean(axis=0)
        std = occupations.std(axis=0, ddof=1)
        noise_of_mean = std / np.sqrt(len(poisson_data))
        signal = qm - INITIAL
        error = mean - qm

        qa, qr, qc = quality(signal, error, noise_of_mean)
        intervals = {}
        for label, start, end in (
            ("0_50", 0.0, 50.0),
            ("50_100", 50.0, 100.0),
            ("100_150", 100.0, 150.0),
            ("0_150", 0.0, 150.0),
            ("150_200", 150.0, 200.0),
            ("0_200", 0.0, 200.0),
        ):
            wqa, wqr, wqc = interval_quality(times, signal, error, noise_of_mean, start, end)
            intervals[f"Q_accuracy_{label}"] = wqa
            intervals[f"Q_repeat_{label}"] = wqr
            intervals[f"Q_{label}"] = wqc

        min_q_to_150 = math.inf
        for start in range(0, 150, 10):
            end = start + 10
            wqa, wqr, wqc = interval_quality(times, signal, error, noise_of_mean, start, end)
            min_q_to_150 = min(min_q_to_150, wqc)
            window_rows.append(
                {
                    "stage": stage,
                    "case_id": case_id,
                    "window": f"{start}-{end}",
                    "wc_eV": float(base["wc_eV"]),
                    "eta": float(base["eta"]),
                    "back_replicas": int(base["back_replicas"]),
                    "Q_accuracy": wqa,
                    "Q_repeat": wqr,
                    "Q_conservative": wqc,
                }
            )

        ntraj = int(base["ntraj"])
        repeats = len(case_rows)
        q0150 = float(intervals["Q_0_150"])
        metrics.append(
            {
                "stage": stage,
                "case_id": case_id,
                "param_id": param_id,
                "origin": base["origin"],
                "random_seed": int(base["random_seed"]),
                "wc_eV": float(base["wc_eV"]),
                "eta": float(base["eta"]),
                "delE_eV": float(base["delE_eV"]),
                "ntraj_per_repeat": ntraj,
                "independent_repeats": repeats,
                "total_forward_paths": ntraj * repeats,
                "back_replicas": int(base["back_replicas"]),
                "tmax": float(times[-1]),
                "signal_rms_0_200": rms(signal),
                "poisson_rms_vs_qm_0_200": rms(error),
                "repeat_noise_of_mean_0_200": rms(noise_of_mean),
                "Q_accuracy_0_200_all": qa,
                "Q_repeat_0_200_all": qr,
                "Q_0_200_all": qc,
                **intervals,
                "min_10au_window_Q_0_150": min_q_to_150,
                "passes_Q0_150_gt_10": q0150 > 10.0,
                "projected_ntraj_per_repeat_Q10_0_150": int(
                    math.ceil(ntraj * (10.0 / max(q0150, 1e-30)) ** 2)
                ),
                "occupation_min": float(occupations.min()),
                "occupation_max": float(occupations.max()),
                "max_particle_error": float(np.max(np.abs(occupations.sum(axis=2) - 5.0))),
            }
        )
        cached[case_id] = (times, mean, std, qm)

    if not metrics:
        raise RuntimeError(f"no complete cases for {stage}")
    metrics.sort(key=lambda row: (float(row["Q_0_150"]), float(row["Q_0_200"])), reverse=True)
    for rank, row in enumerate(metrics, 1):
        row["rank"] = rank
    field_order = ["rank"] + [key for key in metrics[0] if key != "rank"]
    metrics = [{key: row[key] for key in field_order} for row in metrics]

    write_rows(SCAN / f"{stage}_metrics.csv", metrics)
    write_rows(SCAN / f"{stage}_window_metrics.csv", window_rows)
    (SCAN / f"{stage}_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    plot_root = SCAN / f"{stage}_plots"
    plot_root.mkdir(exist_ok=True)
    for row in metrics:
        case_id = str(row["case_id"])
        times, mean, std, qm = cached[case_id]
        fig, axes = plt.subplots(2, 5, figsize=(18, 7.6), sharex=True)
        for orbital, axis in enumerate(axes.flat):
            offset = INITIAL[orbital]
            axis.fill_between(
                times,
                mean[:, orbital] - offset - std[:, orbital],
                mean[:, orbital] - offset + std[:, orbital],
                color="#D55E00",
                alpha=0.13,
                label="repeat spread" if orbital == 0 else None,
            )
            axis.plot(times, mean[:, orbital] - offset, color="#D55E00", lw=1.0,
                      label="Poisson mean" if orbital == 0 else None)
            axis.plot(times, qm[:, orbital] - offset, color="#0072B2", lw=1.3,
                      label="QM" if orbital == 0 else None)
            axis.axvline(150.0, color="#6A3D9A", lw=0.8, ls="--")
            axis.axhline(0.0, color="#aaaaaa", lw=0.55)
            axis.grid(color="#e3e3e3", lw=0.45)
            axis.set_title(f"orbital {orbital}")
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        axes[0, 0].legend(loc="best", fontsize=8)
        for axis in axes[-1]:
            axis.set_xlabel("time (a.u.)")
        for axis in axes[:, 0]:
            axis.set_ylabel(r"$\Delta n_i(t)$")
        fig.suptitle(
            f"rank {row['rank']} | wc={row['wc_eV']:g} eV, eta={row['eta']:.7g}, "
            f"Q(0-150)={row['Q_0_150']:.3g}\n"
            f"Ntraj/run={int(row['ntraj_per_repeat']):,}; repeats={row['independent_repeats']}; "
            f"total forward={int(row['total_forward_paths']):,}; "
            f"back replicas={row['back_replicas']}",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.925))
        fig.savefig(plot_root / f"rank{int(row['rank']):02d}_{case_id}.png", dpi=175)
        plt.close(fig)

    controls = [row for row in metrics if row["origin"] == "replica_control"]
    if controls:
        controls.sort(key=lambda row: int(row["back_replicas"]))
        fig, axis = plt.subplots(figsize=(8.8, 5.4))
        x = [int(row["back_replicas"]) for row in controls]
        axis.plot(x, [float(row["Q_0_150"]) for row in controls], "o-", label="Q(0-150)")
        axis.plot(x, [float(row["Q_100_150"]) for row in controls], "s-", label="Q(100-150)")
        axis.axhline(10.0, color="#6A3D9A", ls="--", lw=1.0, label="target Q=10")
        axis.set_xlabel("backward replicas")
        axis.set_ylabel("conservative Q")
        axis.set_title(
            "Backward-replica pilot | wc=4 eV, eta=1.7e-4\n"
            "Ntraj/run=8,000,000; repeats=3; total forward=24,000,000"
        )
        axis.grid(color="#e3e3e3", lw=0.6)
        axis.legend()
        fig.tight_layout()
        fig.savefig(SCAN / f"{stage}_back_replica_effect.png", dpi=190)
        plt.close(fig)

    summary = {
        "stage": stage,
        "cases": len(metrics),
        "target": "Q(0-150)>10 with tmax=200 and ntraj/run<=24,000,000",
        "best": metrics[0],
    }
    (SCAN / f"{stage}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    analyze(sys.argv[1] if len(sys.argv) > 1 else "pilot")
