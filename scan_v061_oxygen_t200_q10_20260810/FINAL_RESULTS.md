# v0.61 oxygen t=200 parameter scan

## Target

- Output range: 0-200 a.u.
- Required interval: 0-150 a.u.
- Required conservative quality: Q > 10
- Per-run forward path limit: 24,000,000
- Fixed delE: -17.194874968839816 eV

## Best measured case

- Case: f01_random_promoted (random seed 20260817)
- wc: 3.351907 eV
- eta: 3.1184e-05
- Ntraj/run: 16,000,000
- Repeats: 3
- Total forward paths used for the three-repeat estimate: 48,000,000
- Backward replicas: 32
- Q(0-150): 11.244155
- Q_accuracy(0-150): 11.244155
- Q_repeat(0-150): 11.358408
- Q(100-150): 6.761617
- Q(150-200): 2.987081
- Projected Ntraj/run for Q=10: 12,655,123
- Target achieved: True

## Conclusion

The best measured conservative Q is 11.244155, exceeding the
Q=10 target by 1.244155, while using
16,000,000 forward paths per run. This is below the
24,000,000 path cap. The 100-150 a.u. segment remains the weakest interval
at Q=6.761617; the stated target applies to the aggregate
0-150 a.u. interval rather than every sub-window.

## Files

- `high_path_final_metrics.csv`: formal high-path validation candidates
- `all_stage_metrics_final.csv`: all 74 evaluated cases
- `high_path_final_comparison.png`: high-path Q comparison
- `confirm4_plots/rank01_f01_random_promoted.png`: all-orbital best-case plot
- `high_path_runtime.csv`: actual scheduler runtimes
