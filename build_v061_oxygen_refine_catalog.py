from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT / "scan_v061_oxygen_grid_20260809" / "oxygen_grid_combination_catalog.csv"
SCAN = ROOT / "scan_v061_oxygen_local_refine_20260810"
OUTPUT = SCAN / "refinement_combination_catalog.csv"

FIELDS = [
    "catalog_id",
    "generation",
    "selection_status",
    "version",
    "case_id",
    "wc_eV",
    "eta",
    "delE_eV",
    "ntraj_per_repeat",
    "independent_repeats",
    "tmax",
    "qm_signal_rms",
    "poisson_rms_vs_qm",
    "repeat_noise_of_mean_rms",
    "Q_accuracy",
    "Q_repeat",
    "Q_conservative",
    "Q_0_40",
    "Q_40_end",
    "Q_0_60",
    "min_window_Q",
    "projected_ntraj_Q10",
    "selection_score",
    "rank",
    "source",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def normalized(row: dict[str, str], generation: str, retained: set[str]) -> dict[str, str]:
    result = {field: row.get(field, "") for field in FIELDS}
    result["catalog_id"] = f"{generation}_{row.get('case_id', row.get('catalog_id', 'unknown'))}"
    result["generation"] = generation
    if generation == "previous_grid":
        result["catalog_id"] = row.get("catalog_id", result["catalog_id"])
        result["selection_status"] = row.get("selection_status", "reference")
        result["source"] = row.get("source", "v0.61 previous coarse/fine grid")
    else:
        result["selection_status"] = (
            "retained_top30" if row.get("case_id", "") in retained else "screened_out"
        )
        result["source"] = (
            "v0.61 local 3x2M screen, t=0-75 a.u."
            if generation == "local_screen"
            else "v0.61 local 3x8M confirmation, t=0-100 a.u."
        )
    return result


def main() -> None:
    rows: list[dict[str, str]] = []
    rows.extend(normalized(row, "previous_grid", set()) for row in read_rows(PREVIOUS))

    for stage, generation in (("screen", "local_screen"), ("confirm", "local_confirm")):
        metrics = read_rows(SCAN / f"{stage}_metrics.csv")
        retained = {
            row["case_id"] for row in read_rows(SCAN / f"{stage}_retained_top30.csv")
        }
        rows.extend(normalized(row, generation, retained) for row in metrics)

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
