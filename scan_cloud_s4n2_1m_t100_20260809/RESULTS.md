# Four-orbital, two-electron convergence scan

## Fixed setup

- Model: Poisson many-body run (`JOB=1`)
- Orbitals/electrons: `Norb=4`, `Nel=2`
- Sampling: 1,000,000 trajectories per independent run
- Time grid: 0 to 100 a.u., `dt=0.5`, 201 measured points
- Electronic offset: `AHM_DELE_EV=-2.5`
- Convergence windows: 10 a.u.
- Target: strict `Q > 10` in every window through 100 a.u.
- Active-orbital rule: at least 5% of the largest signal before 40 a.u.
- Independent repeats: three for exploration and six for validated finalists

The full search contains 33 distinct `(wc, eta)` combinations and 129 independent
runs. The explored range was `0.05 <= wc <= 5 eV` and
`1e-5 <= eta <= 1.2e-4`; the final local refinement focused on
`0.06 <= wc <= 0.11 eV` and `3e-5 <= eta <= 5e-5`.

## Validated result

The best result after six independent repeats is:

- `wc = 0.25 eV`
- `eta = 6e-5`
- strict `Q(0..100) = 2.1841`
- active orbitals: 0, 2, 3
- inactive orbital: 1
- maximum particle-number error: `1.65e-12`

Strict Q by time window:

| Window (a.u.) | Q |
|---:|---:|
| 0-10 | 113.962 |
| 10-20 | 32.590 |
| 20-30 | 13.150 |
| 30-40 | 9.188 |
| 40-50 | 5.429 |
| 50-60 | 4.469 |
| 60-70 | 3.750 |
| 70-80 | 2.911 |
| 80-90 | 2.240 |
| 90-100 | 2.184 |

Thus one million trajectories satisfy `Q > 10` only through 30 a.u. The fitted
late-time behavior is approximately `Q(t) proportional to t^(-1.40)` over the
25-95 a.u. window centers. All leading parameter combinations show the same
long-time falloff and cluster near `Q=2` at 100 a.u.

Assuming ordinary Monte Carlo scaling, `Q` grows as the square root of the
trajectory count. The best result therefore requires approximately

`1,000,000 * (10 / 2.1841)^2 = 20,962,684`

trajectories to reach `Q=10` through 100 a.u. This is an extrapolation and should
be verified with a larger run before using it as a production guarantee.

## Main conclusion

Changing `wc` and `eta` changes the physical signal amplitude and the ranking in
individual time windows, but it did not remove the accumulated stochastic
variance. A three-repeat apparent winner at `wc=0.11 eV`, `eta=3.5e-5` fell from
`Q=2.47` to `Q=2.09` after six repeats. The stable result is therefore a broad
parameter plateau rather than a sharp optimum.

Further refinement of only `wc` and `eta` is unlikely to produce `Q > 10` through
100 a.u. at one million trajectories. Reaching that target requires about 21
million effective trajectories at the present estimator, or a variance-reduction
change in the Poisson sampling/weight estimator.

## Reproducible outputs

- Final ranking: `round4/strict_final_q100.csv`
- Compact trend table: `round4/trend_summary.csv`
- All-run manifest: `round4/manifest_final.csv`
- Q history: `round4/plots/top_q_by_window.png`
- Best per-orbital result: `round4/plots/best_all_orbitals.png`
