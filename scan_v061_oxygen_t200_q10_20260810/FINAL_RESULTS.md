# v0.61 oxygen t=200 parameter scan

## Target

- Output range: 0-200 a.u.
- Required interval: 0-150 a.u.
- Required conservative quality: Q > 10
- Per-run forward path limit: 24,000,000
- Fixed delE: -17.194874968839816 eV

## Best measured case

- Case: d01_random (random seed 20260816)
- wc: 3.389460 eV
- eta: 2.961e-05
- Ntraj/run: 24,000,000
- Repeats: 3
- Total forward paths used for the three-repeat estimate: 72,000,000
- Backward replicas: 12
- Q(0-150): 9.400467
- Q_accuracy(0-150): 9.868264
- Q_repeat(0-150): 9.400467
- Q(100-150): 5.930639
- Q(150-200): 2.778092
- Projected Ntraj/run for Q=10: 27,158,915
- Target achieved: False

## Conclusion

The best measured conservative Q is 9.400467, which is below 10 by
0.599533. Parameter tuning reduced the projected
per-run requirement to 27,158,915, but
that still exceeds the 24,000,000 path cap. The remaining limitation is the
long-time forward-path variance, especially in 100-150 a.u.; further random
wc/eta tuning is not justified without another variance-reduction change.

## Files

- `high_path_final_metrics.csv`: six 24M-path candidates
- `all_stage_metrics_final.csv`: all 68 evaluated cases
- `high_path_final_comparison.png`: high-path Q comparison
- `confirm2_plots/rank01_d01_random.png`: all-orbital best-case plot
- `high_path_runtime.csv`: actual scheduler runtimes
