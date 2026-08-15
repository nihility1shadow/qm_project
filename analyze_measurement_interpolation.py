#!/usr/bin/env python3
"""Quantify interpolation error from a SepMB measurement schedule on QM data."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from analyze_short_qm_accuracy import load_dat


ACTIVE = np.asarray((0, 5, 6, 7, 8, 9), dtype=int)


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def measurement_steps(nwf: int, period_steps: int, forced_stride: int) -> list[int]:
    dense_end = min(nwf, max(16, period_steps // 8))
    mid_end = min(nwf, max(dense_end, period_steps // 4))
    slow_end = min(nwf, max(mid_end, period_steps // 2))
    stride_mid = max(1, period_steps // 128)
    stride_slow = max(stride_mid, period_steps // 64)
    stride_late = max(stride_slow, period_steps // 32)
    steps = [0]
    j = 1
    while j <= nwf:
        stride = forced_stride if forced_stride > 0 else (
            1 if j <= dense_end else stride_mid if j <= mid_end
            else stride_slow if j <= slow_end else stride_late
        )
        steps.append(j)
        j += stride
    if steps[-1] != nwf:
        steps.append(nwf)
    return steps


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm", required=True, type=Path)
    parser.add_argument("--period-steps", required=True, type=int)
    parser.add_argument("--stride", required=True, type=int)
    parser.add_argument("--window", default=50.0, type=float)
    parser.add_argument("--target-stop", default=400.0, type=float)
    parser.add_argument("--norb", default=10, type=int)
    args = parser.parse_args()

    qm = load_dat(args.qm)
    times = qm[:, 0]
    dt = float(np.median(np.diff(times)))
    nwf = len(times) - 1
    steps = measurement_steps(nwf, args.period_steps, args.stride)
    occupations = qm[:, 4 : 4 + args.norb]
    interpolated = np.column_stack([
        np.interp(times, times[steps], occupations[steps, orbital])
        for orbital in range(args.norb)
    ])
    signal = occupations - occupations[0]
    error = interpolated - occupations

    print(f"stride={args.stride} nmeas={len(steps)-1} max_step={max(np.diff(steps))}")
    start = 0.0
    minimum_q = math.inf
    while start < args.target_stop - 1.0e-12:
        stop = min(start + args.window, args.target_stop)
        mask = (times >= start - 1.0e-12) & (times <= stop + 1.0e-12)
        error_rms = rms(error[mask][:, ACTIVE])
        q = rms(signal[mask][:, ACTIVE]) / error_rms if error_rms else math.inf
        minimum_q = min(minimum_q, q)
        print(f"{start:g}-{stop:g}: interpolation_Q={q:.7g} error_rms={error_rms:.7e}")
        start = stop
    print(f"minimum_interpolation_Q_0_{args.target_stop:g}={minimum_q:.7g}")
    print(f"max_absolute_interpolation_error={np.max(np.abs(error)):.7e}")


if __name__ == "__main__":
    main()
