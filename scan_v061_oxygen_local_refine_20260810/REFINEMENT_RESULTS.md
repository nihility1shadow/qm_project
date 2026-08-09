# v0.61 oxygen 2p local parameter refinement

## Fixed physical configuration

- `Norb = 10`, `Nel = 5`
- `delE = -17.194874968839816 eV = -0.6319 Ha`
- `delE` is the fixed oxygen atomic outermost `2p` orbital energy.
- Only `wc` and `eta` were scanned. The source and executable were not changed.

## Local screen

The screen used:

- `wc = 2, 3, 4, 5, 6, 8 eV`
- `eta = 1.4e-4, 1.7e-4, 2.0e-4, 2.3e-4, 2.6e-4`
- 30 combinations
- `3 x 2,000,000` trajectories per combination
- `t = 0...75 a.u.`, `dt = 0.5 a.u.`
- 32 MPI ranks per Poisson run

The previous winner, `wc=4 eV, eta=2.0e-4`, remained rank 1:

| metric | value |
|---|---:|
| total conservative Q | 4.581 |
| Q, 0-40 a.u. | 13.347 |
| Q, 40-75 a.u. | 3.035 |
| Q, 0-60 a.u. | 7.761 |
| projected trajectories for total Q=10 | 9,530,694 |

This confirms that the previous one-million-path result was not an isolated favorable seed.

## High-path confirmation

The top 30% from the screen were rerun with:

- 9 retained combinations
- `3 x 8,000,000` trajectories per combination
- `t = 0...100 a.u.`, `dt = 0.5 a.u.`
- 64 MPI ranks per Poisson run

The best long-time combination shifted to a lower `eta`:

| rank | wc (eV) | eta | total Q | Q 0-40 | Q 40-100 | Q 0-60 | projected paths for Q=10 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4 | 1.7e-4 | 4.194 | 24.198 | 3.129 | 13.074 | 45,484,672 |
| 2 | 5 | 1.7e-4 | 4.171 | 25.089 | 3.126 | 14.275 | 45,989,891 |
| 3 | 6 | 2.3e-4 | 3.546 | 26.087 | 2.664 | 12.638 | 63,622,032 |
| 4 | 5 | 2.0e-4 | 3.702 | 25.630 | 2.771 | 12.853 | 58,382,894 |
| 7 | 4 | 2.0e-4 | 3.273 | 25.738 | 2.439 | 11.363 | 74,693,877 |

The rank-1 and rank-2 results are statistically close. The defensible parameter region is therefore

`wc = 4...5 eV`, `eta approximately 1.7e-4`,

not a uniquely resolved value of `wc`.

## Interpretation

1. `eta` controls the effective hybridization and the Monte Carlo weight variance more strongly than `wc` in this region.
2. Increasing `eta` strengthens the exact-QM occupation signal, but it also increases late-time path-weight dispersion faster than it increases the signal.
3. Lowering `eta` from `2.0e-4` to `1.7e-4` reduces the exact-QM RMS signal by about 16%, while improving the 60-100 a.u. agreement enough to win the total score.
4. The broad `wc=4...5 eV` plateau is preferable to selecting a single noisy grid point.
5. The result is well converged through roughly 60 a.u. At later times the mean still follows the QM scale, but the repeat band expands and limits the total Q.

The projected `45.5 million` trajectories is a same-time (`t=100 a.u.`) Monte Carlo estimate. It must not be extrapolated directly to `t=4000 a.u.`, because path-weight variance grows with propagation time.

## Runtime

| stage | trajectories/run | time range | MPI ranks | measured run time |
|---|---:|---:|---:|---:|
| local screen | 2,000,000 | 0-75 a.u. | 32 | mean 69.7 s, range 41-140 s |
| confirmation | 8,000,000 | 0-100 a.u. | 64 | mean 142.7 s, range 124-169 s |

Queue waiting time is excluded.

## Files

- `refinement_combination_catalog.csv`: previous 66 rows plus 30 local-screen and 9 confirmation rows.
- `screen_metrics.csv`, `confirm_metrics.csv`: combination-level metrics.
- `screen_window_metrics.csv`, `confirm_window_metrics.csv`: 10-a.u. window metrics.
- `screen_parameter_maps.png`, `confirm_parameter_maps.png`: parameter maps.
- `screen_top_plots/`, `confirm_top_plots/`: all-orbital plots for retained candidates.
