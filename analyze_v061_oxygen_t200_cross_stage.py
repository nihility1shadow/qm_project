from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SCAN = ROOT / "scan_v061_oxygen_t200_q10_20260810"
INITIAL = np.r_[np.ones(5), np.zeros(5)]


def load(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.shape[1] != 14 or not np.isfinite(data).all():
        raise ValueError(f"invalid data: {path}")
    return data


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def quality(signal: np.ndarray, error: np.ndarray, noise: np.ndarray) -> dict[str, float]:
    signal_rms = rms(signal)
    qa = signal_rms / max(rms(error), 1e-30)
    qr = signal_rms / max(rms(noise), 1e-30)
    return {"Q_accuracy": qa, "Q_repeat": qr, "Q_conservative": min(qa, qr)}


def main() -> None:
    all_metrics: list[dict[str, str]] = []
    for stage in ("pilot", "screen", "refine", "fine"):
        path = SCAN / f"{stage}_metrics.csv"
        with path.open(encoding="utf-8-sig") as handle:
            all_metrics.extend(csv.DictReader(handle))
    all_metrics.sort(key=lambda row: float(row["Q_0_150"]), reverse=True)
    fieldnames = list(all_metrics[0])
    with (SCAN / "screening_through_fine_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_metrics)

    repeated = (
        ("refine", "r02_grid"),
        ("fine", "f05_grid"),
    )
    poisson: list[np.ndarray] = []
    batch_q: list[dict[str, object]] = []
    qm_data: np.ndarray | None = None
    for stage, case_id in repeated:
        submission_path = SCAN / f"{stage}_submission.csv"
        with submission_path.open(encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        case_rows = sorted(
            (row for row in rows if row["case_id"] == case_id and row["method"] == "v061"),
            key=lambda row: int(row["repeat"]),
        )
        for row in case_rows:
            poisson.append(
                load(
                    SCAN
                    / f"{stage}_data"
                    / case_id
                    / "v061"
                    / f"rep{row['repeat']}"
                    / f"ahm-sepmb-s10-n5-{row['ntraj']}.dat"
                )
            )
        current_qm = load(SCAN / f"{stage}_data" / "qm" / case_id / "ahm-qm-s10-n5.dat")
        if qm_data is None:
            qm_data = current_qm
        elif not np.allclose(qm_data, current_qm, rtol=0.0, atol=2e-12):
            raise ValueError("matching-parameter QM files differ")
        metric_path = SCAN / f"{stage}_metrics.csv"
        with metric_path.open(encoding="utf-8-sig") as handle:
            row = next(item for item in csv.DictReader(handle) if item["case_id"] == case_id)
        batch_q.append({"stage": stage, "case_id": case_id, "Q_0_150": float(row["Q_0_150"])})

    assert qm_data is not None
    times = qm_data[:, 0]
    occupations = np.stack([data[:, 4:14] for data in poisson])
    qm = qm_data[:, 4:14]
    mean = occupations.mean(axis=0)
    std = occupations.std(axis=0, ddof=1)
    noise = std / np.sqrt(len(poisson))
    signal = qm - INITIAL
    error = mean - qm

    windows: dict[str, dict[str, float]] = {}
    for label, start, end in (
        ("0_50", 0.0, 50.0),
        ("50_100", 50.0, 100.0),
        ("100_150", 100.0, 150.0),
        ("0_150", 0.0, 150.0),
        ("150_200", 150.0, 200.0),
        ("0_200", 0.0, 200.0),
    ):
        mask = (times >= start) & (times <= end)
        windows[label] = quality(signal[mask], error[mask], noise[mask])

    result = {
        "wc_eV": 3.4,
        "eta": 4e-5,
        "delE_eV": -17.194874968839816,
        "ntraj_per_repeat": 2_000_000,
        "independent_repeats": len(poisson),
        "total_forward_paths": len(poisson) * 2_000_000,
        "back_replicas": 12,
        "stratify_forward": 1,
        "batch_Q_0_150": batch_q,
        "combined_windows": windows,
        "projected_ntraj_per_repeat_Q10_0_150": int(
            np.ceil(2_000_000 * (10.0 / windows["0_150"]["Q_conservative"]) ** 2)
        ),
        "max_particle_error": float(np.max(np.abs(occupations.sum(axis=2) - 5.0))),
    }
    (SCAN / "cross_stage_replicate_check.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(2, 5, figsize=(18, 7.6), sharex=True)
    for orbital, axis in enumerate(axes.flat):
        offset = INITIAL[orbital]
        axis.fill_between(
            times,
            mean[:, orbital] - offset - std[:, orbital],
            mean[:, orbital] - offset + std[:, orbital],
            color="#D55E00",
            alpha=0.13,
        )
        axis.plot(times, mean[:, orbital] - offset, color="#D55E00", lw=1.0)
        axis.plot(times, qm[:, orbital] - offset, color="#0072B2", lw=1.3)
        axis.axvline(150.0, color="#6A3D9A", lw=0.8, ls="--")
        axis.axhline(0.0, color="#aaaaaa", lw=0.55)
        axis.grid(color="#e3e3e3", lw=0.45)
        axis.set_title(f"orbital {orbital}")
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    for axis in axes[-1]:
        axis.set_xlabel("time (a.u.)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\Delta n_i(t)$")
    fig.suptitle(
        f"Cross-stage repeat check | wc=3.4 eV, eta=4e-5, "
        f"Q(0-150)={windows['0_150']['Q_conservative']:.3g}\n"
        f"Ntraj/run=2,000,000; repeats={len(poisson)}; "
        f"total forward={len(poisson) * 2_000_000:,}; back replicas=12",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    fig.savefig(SCAN / "cross_stage_wc3p4_eta4e-5_six_repeats.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
