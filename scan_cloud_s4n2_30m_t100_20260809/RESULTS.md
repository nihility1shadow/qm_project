# Two-run 30-million-trajectory result

## Run configuration

- Orbitals/electrons: `Norb=4`, `Nel=2`
- Trajectories: 30,000,000 per independent run
- Independent runs: 2
- Parameters: `wc=0.25 eV`, `eta=6e-5`, `AHM_DELE_EV=-2.5`
- Time range: 0 to 100 a.u., 201 measured points at 0.5 a.u. intervals
- MPI ranks: 64 per run
- Jobs: `636247`, `636248`
- Elapsed time: 3:57 and 3:50

## Strict convergence result

Using the same 10 a.u. windows and active-orbital definition as the one-million
trajectory scan, the strict single-run convergence estimate is:

`Q(0..100) = 9.4902`

| Window (a.u.) | Q |
|---:|---:|
| 0-10 | 614.139 |
| 10-20 | 206.311 |
| 20-30 | 80.895 |
| 30-40 | 50.003 |
| 40-50 | 30.984 |
| 50-60 | 25.434 |
| 60-70 | 18.536 |
| 70-80 | 15.575 |
| 80-90 | 12.401 |
| 90-100 | 9.490 |

The bottleneck is orbital 3 in the 90-100 a.u. window. Orbital Q values in that
window are 15.601, 17.599, 16.266, and 9.490 for orbitals 0 through 3.

The maximum particle-number error is `5.09e-11`, so the result remains consistent
with two-electron number conservation.

## Scaling interpretation

The validated one-million-trajectory baseline was `Q=2.1841`. Increasing the
trajectory count by a factor of 30 improved the final strict Q by a factor of
4.345, compared with the ideal square-root prediction `sqrt(30)=5.477`.

Most windows before 90 a.u. are close to square-root scaling. The shortfall is
concentrated in orbital 3 during the final window and is particularly uncertain
because only two independent repeats are available.

At the observed scaling, a single run would need approximately

`30,000,000 * (10 / 9.4902)^2 = 33,309,754`

trajectories to reach strict `Q=10` through 100 a.u. If the two 30-million runs
are averaged and noise is interpreted as the standard error of that mean, their
effective 60-million-trajectory estimate is `Q_mean=13.42`. This second value
passes 10, but its uncertainty estimate has only one degree of freedom; a third
independent run is needed for a robust convergence claim.

## Files

- Raw-run manifest: `manifest.csv`
- Strict ranking: `strict_30m_two_repeats_q100.csv`
- Window scaling table: `plots/q_scaling_by_window.csv`
- Q scaling plot: `plots/q_1m_vs_30m.png`
- Per-orbital plot: `plots/best_all_orbitals.png`
