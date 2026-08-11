#!/usr/bin/env python3
"""Combine independent SepMB data files with trajectory-count weights."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_run(value: str) -> tuple[int, Path]:
    fields = value.split("=", 1)
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("run must use NTRAJ=PATH")
    return int(fields[0]), Path(fields[1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, action="append", type=parse_run)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    arrays: list[np.ndarray] = []
    weights: list[int] = []
    headers: list[str] = []
    for index, (ntraj, path) in enumerate(args.run):
        data = np.loadtxt(path, comments="#")
        if data.ndim == 1:
            data = data[None, :]
        if arrays:
            if data.shape != arrays[0].shape or not np.allclose(data[:, 0], arrays[0][:, 0]):
                raise ValueError(f"incompatible time grid: {path}")
        else:
            with path.open(encoding="utf-8") as stream:
                headers = [line.rstrip() for line in stream if line.startswith("#")]
        arrays.append(data)
        weights.append(ntraj)

    total = sum(weights)
    combined = np.zeros_like(arrays[0])
    combined[:, 0] = arrays[0][:, 0]
    for data, weight in zip(arrays, weights):
        combined[:, 1:] += data[:, 1:] * (weight / total)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"#COMBINED_INDEPENDENT_RUNS: count={len(arrays)} total_ntraj={total}\n")
        for line in headers:
            stream.write(line + "\n")
        np.savetxt(
            stream,
            combined,
            fmt=["%12.8f"] + ["%+1.16e"] * (combined.shape[1] - 1),
        )
    print(f"{args.output}: {len(arrays)} runs, total_ntraj={total}")


if __name__ == "__main__":
    main()
