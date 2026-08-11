#!/usr/bin/env python3
"""Merge v0.71 coarse, refine, confirmation, and final Q tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


STAGES = ("coarse", "refine", "confirm", "final")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    fields: list[str] = ["stage"]
    for stage in STAGES:
        path = args.root / f"{stage}-analysis" / "grid_metrics.csv"
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8-sig") as stream:
            stage_rows = list(csv.DictReader(stream))
        for row in stage_rows:
            rows.append({"stage": stage, **row})
            for field in row:
                if field not in fields:
                    fields.append(field)

    rows.sort(
        key=lambda row: float(row.get("min_Q_active_0_150", "-inf")),
        reverse=True,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"rows={len(rows)} output={args.output}")


if __name__ == "__main__":
    main()
