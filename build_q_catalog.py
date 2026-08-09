#!/usr/bin/env python3
"""Combine historical and current strict-Q rankings into one lookup catalog."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


HISTORICAL_SOURCES = [
    ("scan_cloud_strict_t50_20260809/strict_per_orbital.csv", "pre-resonant", "legacy-t50"),
    ("scan_cloud_strict_refined_t50_20260809/strict_per_orbital.csv", "pre-resonant", "legacy-refined-t50"),
    ("scan_cloud_strict_validation_t75_20260809/strict_per_orbital_q60.csv", "pre-resonant", "legacy-validation-t75"),
    ("scan_cloud_dense_ab_20260809/strict_ab_q60.csv", "original-delE", "dense-measure-ab"),
    ("scan_cloud_grid_1m_q4_20260809/strict_grid_q60.csv", "original-delE", "original-coarse-1m"),
    ("scan_cloud_grid_1m_fine_20260809/strict_fine_q60.csv", "original-delE", "original-fine-1m"),
    ("scan_cloud_resonant_probe_1m_20260809/strict_resonant_q60.csv", "delE=-2.5eV", "resonant-probe-1m"),
    ("scan_cloud_resonant_fine_1m_20260809/strict_resonant_fine_q60.csv", "delE=-2.5eV", "resonant-fine-1m"),
    ("scan_cloud_resonant_micro_1m_20260809/strict_resonant_micro_grid3_q60.csv", "delE=-2.5eV", "resonant-micro-grid3-1m"),
    ("scan_cloud_resonant_micro_1m_20260809/strict_resonant_micro_validated8_q60.csv", "delE=-2.5eV", "resonant-micro-validated8-1m"),
    ("scan_cloud_resonant_lowwc_1m_20260809/strict_resonant_lowwc_q60.csv", "delE=-2.5eV", "resonant-lowwc-1m"),
    ("scan_cloud_resonant_verylowwc_1m_20260809/strict_resonant_verylowwc_q60.csv", "delE=-2.5eV", "resonant-verylowwc-1m"),
    ("scan_cloud_resonant_ultralowwc_1m_20260809/strict_resonant_ultralowwc_grid3_q60.csv", "delE=-2.5eV", "resonant-ultralowwc-grid3-1m"),
    ("scan_cloud_resonant_ultralowwc_1m_20260809/strict_resonant_ultralowwc_validated8_q60.csv", "delE=-2.5eV", "resonant-ultralowwc-validated8-1m"),
]

WINDOWS = [(start, start + 10) for start in range(0, 150, 10)]


def number(value: str | None) -> float | int | str:
    if value is None or value == "":
        return ""
    try:
        parsed = float(value)
    except ValueError:
        return value
    if parsed.is_integer() and "e" not in value.lower() and "." not in value:
        return int(parsed)
    return parsed


def band_class(wc: float) -> str:
    if wc < 1.0:
        return "ultranarrow"
    if wc < 4.0:
        return "narrow"
    return "wide"


def selection_metadata(ranking_path: Path) -> dict[str, dict[str, str]]:
    selection_path = ranking_path.parent / "selection.csv"
    if not selection_path.exists():
        return {}
    with selection_path.open(newline="", encoding="utf-8-sig") as handle:
        return {row["case_id"]: row for row in csv.DictReader(handle)}


def source_rows(root: Path, path: Path, branch: str, stage: str) -> list[dict]:
    absolute = root / path
    if not absolute.exists():
        return []
    with absolute.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    selected = selection_metadata(absolute)
    round_match = re.fullmatch(r"round(\d+)", absolute.parent.name)
    search_round = int(round_match.group(1)) if round_match else ""
    output = []
    for row in rows:
        wc = float(row["wc_eV"])
        decision = selected.get(row["case_id"], {})
        record = {
            "source_stage": stage,
            "model_branch": branch,
            "source_file": str(path).replace("\\", "/"),
            "search_round": search_round,
            "selection_rank": number(decision.get("rank")),
            "selection_status": decision.get("selection_status", "not-applicable"),
            "case_id": row["case_id"],
            "wc_eV": wc,
            "eta": float(row["eta"]),
            "band_class": band_class(wc),
            "ntraj_per_run": number(row.get("ntraj_per_run")),
            "independent_repeats": number(row.get("independent_repeats")),
            "time_end": number(row.get("time_end")),
            "required_until": number(row.get("required_until")),
            "required_Q": number(row.get("required_Q")),
            "passes_required_interval": row.get("passes_required_interval", ""),
            "active_orbitals": row.get("active_orbitals", ""),
            "inactive_orbitals": row.get("inactive_orbitals", ""),
            "max_inactive_noise": number(row.get("max_inactive_noise")),
            "inactive_noise_limit": number(row.get("inactive_noise_limit")),
            "signal_rms": number(row.get("signal_rms")),
            "signal_peak": number(row.get("signal_peak")),
            "repeat_noise_rms": number(row.get("repeat_noise_rms")),
            "particle_number_max_error": number(row.get("particle_number_max_error")),
        }
        available = []
        for start, end in WINDOWS:
            key = f"min_active_Q_{start}_{end}"
            value = number(row.get(key))
            record[f"Q_{start}_{end}"] = value
            if isinstance(value, (int, float)) and end <= float(record["required_until"]):
                available.append((float(value), f"{start}-{end}"))
        if available:
            _, record["bottleneck_window"] = min(available)
        else:
            record["bottleneck_window"] = ""
        output.append(record)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()

    sources = [(Path(path), branch, stage) for path, branch, stage in HISTORICAL_SOURCES]
    current = root / "scan_cloud_8m_t150_20260809"
    for ranking in sorted(current.glob("round*/strict*.csv")):
        relative = ranking.relative_to(root)
        stage = ranking.parent.name + "-8m-t150"
        sources.append((relative, "delE=-2.5eV", stage))

    records = []
    for path, branch, stage in sources:
        records.extend(source_rows(root, path, branch, stage))

    fields = [
        "source_stage", "model_branch", "source_file", "search_round",
        "selection_rank", "selection_status", "case_id", "wc_eV",
        "eta", "band_class", "ntraj_per_run", "independent_repeats",
        "time_end", "required_until", "required_Q", "bottleneck_window",
        "passes_required_interval", "active_orbitals", "inactive_orbitals",
        "max_inactive_noise", "inactive_noise_limit", "signal_rms", "signal_peak",
        "repeat_noise_rms", "particle_number_max_error",
    ] + [f"Q_{start}_{end}" for start, end in WINDOWS]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    print(f"rows={len(records)} output={args.output}")


if __name__ == "__main__":
    main()
