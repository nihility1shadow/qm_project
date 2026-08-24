# 2500 a.u. pure-Poisson single-run result

## Constraints

- Algorithm: stochastic Poisson only (`SEP_MB_DETERMINISTIC_FOCK=0`).
- The v0.93 low-memory implementation and all rollback files were preserved.
- No source file was changed on the server during this parameter scan.
- The final result is one MPI run, not an average of independent completed runs.
- System: 10 orbitals, 5 electrons, oxygen outer-orbital energy parameter.

## Selected parameters

| Parameter | Value |
| --- | ---: |
| `AHM_WC_EV` | 2.7 eV |
| `AHM_ETA` | 0.30e-6 |
| `AHM_DELE_EV` | -17.194874968839816 eV |
| `AHM_NSTEP` | 5000 |
| `SEP_MB_TMAX` | 2500 a.u. |
| Forward trajectories | 50,000,000 |
| Back replicas | 256 to 512, linear in time |
| MPI ranks | 384 on 6 nodes |
| Seed | 2500824 |

The narrow eta scan used 0.25e-6, 0.30e-6, and 0.35e-6 with matching QM references. Eta = 0.30e-6 gave the best late-time balance. A separate back-replica scan showed that 512 maximum replicas were worthwhile at this eta, while 1024 replicas cost too much for the gain.

## Final run

- Slurm job: `646440`
- State: completed, exit code 0
- Wall time: 14:45:20
- CPU budget: about 5,666 core-hours
- Output points: 5,001, from 0 to 2500 a.u.
- Peak resident memory: 12,944 KiB per rank, about 4.74 GiB summed over 384 ranks
- Maximum particle-number error: 2.20e-8
- Active-orbital signal RMS: 1.96e-7
- Active-orbital error RMS: 6.81e-8

## Accuracy against matching QM

Q is defined as QM occupation-change RMS divided by Poisson error RMS. Strict Q is the minimum Q among active orbitals 0, 5, 6, 7, 8, and 9 in each 100 a.u. window.

- Minimum aggregate Q over all windows: 1.3046 (2200-2300 a.u.).
- Minimum strict Q over all windows: 0.6949 (2400-2500 a.u., orbital 9).
- Every 100 a.u. window from 0 through 700 a.u. has strict Q >= 10.
- Strict Q remains above 1 through 2200 a.u.; it is not monotonic afterward because the physical signal is very small and the residual Monte Carlo phase noise is random.
- The initially occupied inactive orbitals 1-4 have QM signal RMS only 2.0e-11 to 5.8e-11. Their Poisson RMS absolute error is 7.1e-9 to 8.5e-9, so relative Q is not physically useful for those nearly stationary orbitals.

Full-time active-orbital Q values are 3.00 (orbital 0), 2.54 (5), 2.50 (6), 2.68 (7), 2.21 (8), and 2.21 (9).

## Long-time limit

The late-time limit is the Poisson phase/sign problem, not memory growth. In the 2400-2500 a.u. window, strict Q is 0.6949. Assuming the normal Monte Carlo square-root law:

- Strict Q = 1 would require about 1.04e8 trajectories, or about 31 hours at the measured 384-rank throughput.
- Strict Q = 10 would require about 1.04e10 trajectories, or about 127 days at the same allocation.

Therefore a pure-Poisson, single-run target of strict Q >= 10 over the full 2500 a.u. range cannot be reached within a one-billion-trajectory budget by parameter tuning alone. The low-memory property is retained, but the variance grows as the physical signal becomes smaller than the residual complex-phase cancellation noise.

## Files

- Final data: `final_50m_job_646440/ahm-sepmb-s10-n5-50000000.dat`
- Window metrics: `final_50m_job_646440/analysis/q_windows.csv`
- Full metrics: `final_50m_job_646440/analysis/q_metrics.json`
- All-orbital comparison: `final_50m_job_646440/plots/all_orbitals_poisson_vs_qm.png`
- Q-by-window plot: `final_50m_job_646440/plots/q_by_window.png`
- All screened combinations: `parameter_screen.csv`
