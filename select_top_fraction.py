#!/usr/bin/env python3
"""Label a strict-Q ranking with retained/eliminated status."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--round", required=True, type=int)
    parser.add_argument("--fraction", type=float, default=0.30)
    args = parser.parse_args()

    with args.ranking.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    keep = math.ceil(len(rows) * args.fraction)

    fields = [
        "round", "rank", "selection_status", "case_id", "wc_eV", "eta",
        "required_until", "required_Q", "bottleneck_window",
        "passes_required_interval", "independent_repeats", "ntraj_per_run",
    ]
    output = []
    for row in rows:
        required_until = float(row["required_until"])
        windows = []
        for key, value in row.items():
            if not key.startswith("min_active_Q_") or not value:
                continue
            start, end = map(int, key.removeprefix("min_active_Q_").split("_"))
            if end <= required_until:
                windows.append((float(value), f"{start}-{end}"))
        bottleneck = min(windows)[1] if windows else ""
        rank = int(row["rank"])
        output.append({
            "round": args.round,
            "rank": rank,
            "selection_status": "retained" if rank <= keep else "eliminated",
            "case_id": row["case_id"],
            "wc_eV": row["wc_eV"],
            "eta": row["eta"],
            "required_until": row["required_until"],
            "required_Q": row["required_Q"],
            "bottleneck_window": bottleneck,
            "passes_required_interval": row["passes_required_interval"],
            "independent_repeats": row["independent_repeats"],
            "ntraj_per_run": row["ntraj_per_run"],
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(rows)} retained={keep} output={args.output}")


if __name__ == "__main__":
    main()
