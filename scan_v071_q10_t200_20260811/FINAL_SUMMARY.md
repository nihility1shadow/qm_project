# SepMB v0.71: Q > 10 through 150 a.u.

## Result

The final independent-repeat mean passes the requested aggregate active-orbital
criterion through 150 a.u.:

- Minimum window Q over 0-50, 50-100, and 100-150 a.u.: **10.7760**
- Q(0-50): 109.8492
- Q(50-100): 28.7399
- Q(100-150): 10.7760
- Q(150-200): 4.3853 (reported for extension planning; not part of the target)
- Active-orbital amplitude ratio: 0.997599
- Active-orbital cosine similarity: 0.992572
- Maximum particle-number error: 1.47e-9

Here Q is

`RMS(QM occupation-change signal) / RMS(SepMB - QM error)`.

The aggregate active set is orbitals 0, 5, 6, 7, 8, and 9. Orbitals 1-4
have nearly zero QM signal, so a relative signal/error Q is ill-conditioned for
them. Their absolute results are retained in the all-orbital table and plot.
The stricter minimum of the individual active-orbital Q values through 150 a.u.
is 5.8270; the requested aggregate window criterion is the one that passes 10.

## Physical and numerical parameters

| Parameter | Final value |
|---|---:|
| Physical orbitals, Norb | 10 |
| Electrons, Nel | 5 |
| Oxygen outer-orbital energy, delE | -17.1948749688 eV |
| Bath width, wc | 3.389460 eV |
| Coupling, eta | 2.9610e-5 |
| Time step, dt | 0.5 a.u. |
| Number of time steps | 400 |
| Output time range | 0-200 a.u. |
| Direct output points | 401 |
| Forward trajectories per repeat | 1,250,000 |
| Independent repeats | 2 |
| Total effective forward trajectories | 2,500,000 |
| Backward replicas per forward trajectory, B | 256 |
| MPI ranks per repeat | 128 |
| Concurrent MPI ranks | 256 |

The width remains inside the required range `0 < wc <= 10 eV`.

## Runtime

| Job | Forward trajectories | MPI ranks | Nodes | Elapsed |
|---|---:|---:|---:|---:|
| 645067 | 1,250,000 | 128 | 2 | 14 min 15 s |
| 645068 | 1,250,000 | 128 | 2 | 13 min 42 s |

The two jobs ran concurrently. Therefore, the observed final-validation wall
time was 14 min 15 s, while their summed allocation was about 59.6 core-hours.
The same-parameter QM reference job took 3 min 14 s.

## Search history

The optimization retained the top candidates while adding new random candidates:

1. Coarse screen: 8 combinations, 25,000 trajectories each.
2. Refined screen: 6 combinations, 200,000 trajectories each.
3. Confirmation: 5 combinations, 800,000 trajectories each.
4. Final validation: two independent 1,250,000-trajectory repeats at the best point.

All 20 evaluated rows are preserved in `all_parameter_q_catalog.csv`.

## Files

- Final weighted data: `final/final_best/poisson/ahm-sepmb-s10-n5-2500000.dat`
- Same-parameter QM: `final/final_best/qm/ahm-qm-s10-n5.dat`
- Main plot: `final-analysis/final_best_all_orbitals_qm_vs_sepmb.png`
- Final metrics: `final-analysis/grid_metrics.csv`
- Per-orbital metrics: `final-analysis/all-orbital-metrics/short_qm_accuracy_per_orbital.csv`
- Complete parameter catalog: `all_parameter_q_catalog.csv`

## Reproduction settings

For each independent repeat, use:

```bash
export AHM_WC_EV=3.389460
export AHM_ETA=2.9610e-5
export AHM_DELE_EV=-17.194874968839816
export AHM_NSTEP=400
export SEP_MB_BACK_REPLICAS=256
export SEP_MB_MEASURE_STRIDE=1
mpirun -np 128 ./na_mpi_cloud.out 10 1 1250000 10 5
```

Independent repeats must use independent random seeds. Combine their numeric
columns with trajectory-count weights.
