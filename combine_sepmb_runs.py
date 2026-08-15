#!/usr/bin/env python3
"""Combine independent SepMB data files with trajectory-count weights."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from analyze_short_qm_accuracy import load_dat


def parse_run(value: str) -> tuple[int, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("run must use NTRAJ=PATH")
    raw_weight, raw_path = value.split("=", 1)
    weight = int(raw_weight)
    if weight <= 0:
        raise argparse.ArgumentTypeError("NTRAJ must be positive")
    return weight, Path(raw_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, action="append", type=parse_run)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    arrays: list[np.ndarray] = []
    weights: list[int] = []
    for weight, path in args.run:
        data = load_dat(path)
        if arrays:
            if data.shape != arrays[0].shape:
                raise ValueError(f"shape mismatch: {path} has {data.shape}")
            if not np.allclose(data[:, 0], arrays[0][:, 0], rtol=0.0, atol=1.0e-12):
                raise ValueError(f"time-grid mismatch: {path}")
        arrays.append(data)
        weights.append(weight)

    combined = np.average(np.stack(arrays), axis=0, weights=np.asarray(weights))
    combined[:, 0] = arrays[0][:, 0]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("#COMBINED: independent SepMB trajectory-count weighted mean\n")
        stream.write(f"#total_ntraj={sum(weights)} independent_runs={len(weights)}\n")
        for (weight, path) in args.run:
            stream.write(f"#source ntraj={weight} path={path.resolve()}\n")
        np.savetxt(stream, combined, fmt="%+1.16e")


if __name__ == "__main__":
    main()
