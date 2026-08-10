from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
SCAN = ROOT / "scan_v061_oxygen_t200_q10_20260810"
STAGES = ("pilot", "screen", "refine", "fine", "validate", "confirm", "confirm2")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    all_rows: list[dict[str, str]] = []
    for stage in STAGES:
        all_rows.extend(read_rows(SCAN / f"{stage}_metrics.csv"))
    all_rows.sort(key=lambda row: float(row["Q_0_150"]), reverse=True)
    write_rows(SCAN / "all_stage_metrics_final.csv", all_rows)

    high_path = [row for row in all_rows if row["stage"] in {"confirm", "confirm2"}]
    high_path.sort(key=lambda row: float(row["Q_0_150"]), reverse=True)
    final_rows: list[dict[str, object]] = []
    for rank, row in enumerate(high_path, 1):
        final_rows.append(
            {
                "final_rank": rank,
                "stage": row["stage"],
                "case_id": row["case_id"],
                "origin": row["origin"],
                "random_seed": int(row["random_seed"]),
                "wc_eV": float(row["wc_eV"]),
                "eta": float(row["eta"]),
                "delE_eV": float(row["delE_eV"]),
                "ntraj_per_repeat": int(row["ntraj_per_repeat"]),
                "independent_repeats": int(row["independent_repeats"]),
                "total_forward_paths": int(row["total_forward_paths"]),
                "back_replicas": int(row["back_replicas"]),
                "Q_0_50": float(row["Q_0_50"]),
                "Q_50_100": float(row["Q_50_100"]),
                "Q_100_150": float(row["Q_100_150"]),
                "Q_0_150": float(row["Q_0_150"]),
                "Q_accuracy_0_150": float(row["Q_accuracy_0_150"]),
                "Q_repeat_0_150": float(row["Q_repeat_0_150"]),
                "Q_150_200": float(row["Q_150_200"]),
                "min_10au_window_Q_0_150": float(row["min_10au_window_Q_0_150"]),
                "projected_ntraj_per_repeat_Q10_0_150": int(
                    row["projected_ntraj_per_repeat_Q10_0_150"]
                ),
                "signal_rms_0_200": float(row["signal_rms_0_200"]),
                "max_particle_error": float(row["max_particle_error"]),
                "passes_Q0_150_gt_10": row["passes_Q0_150_gt_10"].lower() == "true",
            }
        )
    write_rows(SCAN / "high_path_final_metrics.csv", final_rows)

    best = final_rows[0]
    runtime = read_rows(SCAN / "high_path_runtime.csv")
    runtime_seconds = [int(row["elapsed_seconds"]) for row in runtime]
    summary = {
        "target": {
            "tmax": 200.0,
            "score_interval": [0.0, 150.0],
            "required_Q": 10.0,
            "maximum_ntraj_per_run": 24_000_000,
        },
        "fixed_physics": {
            "Norb": 10,
            "Nel": 5,
            "delE_eV": -17.194874968839816,
        },
        "sampling": {
            "version": "v0.61",
            "back_replicas": 12,
            "stratify_forward_count": True,
            "independent_repeats": 3,
        },
        "search": {
            "stages": list(STAGES),
            "evaluated_parameter_cases": len(all_rows),
            "high_path_cases": len(final_rows),
            "random_seeds": [20260810, 20260811, 20260812, 20260813, 20260814, 20260815, 20260816],
        },
        "best": best,
        "target_achieved": bool(best["passes_Q0_150_gt_10"]),
        "Q_shortfall": 10.0 - float(best["Q_0_150"]),
        "runtime_seconds": {
            "runs": len(runtime_seconds),
            "minimum": min(runtime_seconds),
            "maximum": max(runtime_seconds),
            "mean": sum(runtime_seconds) / len(runtime_seconds),
        },
        "best_plot": str(SCAN / "confirm2_plots" / "rank01_d01_random.png"),
    }
    (SCAN / "FINAL_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, (axis_q, axis_windows) = plt.subplots(1, 2, figsize=(14, 5.6))
    scatter = axis_q.scatter(
        [float(row["eta"]) * 1e5 for row in final_rows],
        [float(row["Q_0_150"]) for row in final_rows],
        c=[float(row["wc_eV"]) for row in final_rows],
        cmap="viridis",
        s=85,
        edgecolor="#222222",
        linewidth=0.6,
    )
    label_offsets = {
        "d01_random": (5, 6, "left"),
        "c02_random_prev": (-5, 8, "right"),
        "d02_random": (-6, 8, "right"),
        "c04_random_new": (8, 12, "left"),
        "c01_known": (8, -14, "left"),
        "c03_random_new": (-6, 8, "right"),
    }
    for row in final_rows:
        eta_value = float(row["eta"]) * 1e5
        dx, dy, horizontal_alignment = label_offsets[str(row["case_id"])]
        axis_q.annotate(
            str(row["case_id"]),
            (eta_value, float(row["Q_0_150"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8,
            ha=horizontal_alignment,
        )
    axis_q.axhline(10.0, color="#6A3D9A", ls="--", lw=1.1, label="target Q=10")
    axis_q.set_xlabel(r"eta ($10^{-5}$)")
    axis_q.set_ylabel("Q(0-150)")
    axis_q.set_title("High-path parameter comparison")
    axis_q.grid(color="#e2e2e2", lw=0.6)
    axis_q.margins(x=0.08, y=0.08)
    axis_q.legend()
    colorbar = fig.colorbar(scatter, ax=axis_q)
    colorbar.set_label("wc (eV)")

    top = final_rows[:4]
    labels = [str(row["case_id"]) for row in top]
    x = list(range(len(top)))
    width = 0.24
    for offset, key, label, color in (
        (-width, "Q_0_50", "Q(0-50)", "#0072B2"),
        (0.0, "Q_50_100", "Q(50-100)", "#E69F00"),
        (width, "Q_100_150", "Q(100-150)", "#D55E00"),
    ):
        axis_windows.bar(
            [value + offset for value in x],
            [float(row[key]) for row in top],
            width=width,
            label=label,
            color=color,
        )
    axis_windows.axhline(10.0, color="#6A3D9A", ls="--", lw=1.0)
    axis_windows.set_xticks(x, labels, rotation=18, ha="right")
    axis_windows.set_ylabel("conservative Q")
    axis_windows.set_title("Time-window quality for top high-path cases")
    axis_windows.grid(axis="y", color="#e2e2e2", lw=0.6)
    axis_windows.legend()
    fig.suptitle(
        "v0.61 oxygen scan | Ntraj/run=24,000,000; repeats=3; "
        "total forward/case=72,000,000; back replicas=12"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(SCAN / "high_path_final_comparison.png", dpi=190)
    plt.close(fig)

    report = f"""# v0.61 oxygen t=200 parameter scan

## Target

- Output range: 0-200 a.u.
- Required interval: 0-150 a.u.
- Required conservative quality: Q > 10
- Per-run forward path limit: 24,000,000
- Fixed delE: -17.194874968839816 eV

## Best measured case

- Case: {best['case_id']} (random seed {best['random_seed']})
- wc: {best['wc_eV']:.6f} eV
- eta: {best['eta']:.9g}
- Ntraj/run: {best['ntraj_per_repeat']:,}
- Repeats: {best['independent_repeats']}
- Total forward paths used for the three-repeat estimate: {best['total_forward_paths']:,}
- Backward replicas: {best['back_replicas']}
- Q(0-150): {best['Q_0_150']:.6f}
- Q_accuracy(0-150): {best['Q_accuracy_0_150']:.6f}
- Q_repeat(0-150): {best['Q_repeat_0_150']:.6f}
- Q(100-150): {best['Q_100_150']:.6f}
- Q(150-200): {best['Q_150_200']:.6f}
- Projected Ntraj/run for Q=10: {best['projected_ntraj_per_repeat_Q10_0_150']:,}
- Target achieved: {best['passes_Q0_150_gt_10']}

## Conclusion

The best measured conservative Q is {best['Q_0_150']:.6f}, which is below 10 by
{10.0 - float(best['Q_0_150']):.6f}. Parameter tuning reduced the projected
per-run requirement to {best['projected_ntraj_per_repeat_Q10_0_150']:,}, but
that still exceeds the 24,000,000 path cap. The remaining limitation is the
long-time forward-path variance, especially in 100-150 a.u.; further random
wc/eta tuning is not justified without another variance-reduction change.

## Files

- `high_path_final_metrics.csv`: six 24M-path candidates
- `all_stage_metrics_final.csv`: all {len(all_rows)} evaluated cases
- `high_path_final_comparison.png`: high-path Q comparison
- `confirm2_plots/rank01_d01_random.png`: all-orbital best-case plot
- `high_path_runtime.csv`: actual scheduler runtimes
"""
    (SCAN / "FINAL_RESULTS.md").write_text(report, encoding="ascii")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
