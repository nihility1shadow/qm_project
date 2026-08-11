# SepMB v0.71 direct-grid run through 500 a.u.

## Configuration

| Parameter | Value |
|---|---:|
| Norb / Nel | 10 / 5 |
| wc | 3.389460 eV |
| eta | 2.9610e-5 |
| delE | -17.194874968839816 eV |
| dt / nstep / tmax | 0.5 / 1000 / 500 a.u. |
| Direct output points | 1001 |
| Forward paths | 2 x 1,250,000 = 2,500,000 |
| Backward replicas B | 256 |
| MPI ranks | 2 x 128, run concurrently |

The two SepMB repeats used independent random seeds. Their data columns were
combined with forward-trajectory-count weights. No time interpolation was used.

## Runtime

| Job | Role | Elapsed | MPI ranks | Nodes | Exit code |
|---|---|---:|---:|---:|---:|
| 645072 | QM | 6 min 53 s | 1 | 1 | 0:0 |
| 645073 | SepMB repeat 1 | 1 h 48 min 37 s | 128 | 2 | 0:0 |
| 645074 | SepMB repeat 2 | 1 h 50 min 18 s | 128 | 2 | 0:0 |

Because the two SepMB jobs ran concurrently, the final SepMB wall time was
1 h 50 min 18 s. Their summed allocation was about 467 core-hours. The
corresponding 200-a.u. validation took 14 min 15 s, so increasing tmax by 2.5
increased measured wall time by 7.74. This is an empirical comparison, not a
claim of an exact asymptotic exponent.

## Window results

Q is `RMS(QM active-orbital signal) / RMS(SepMB - QM error)` for active
orbitals 0, 5, 6, 7, 8, and 9.

| Window (a.u.) | Q | Amplitude ratio | Cosine | Projected paths for Q=10 |
|---:|---:|---:|---:|---:|
| 0-50 | 48.364 | 1.015 | 0.99990 | 106,878 |
| 50-100 | 21.997 | 1.031 | 0.99946 | 516,688 |
| 100-150 | 10.686 | 1.035 | 0.99635 | 2,189,218 |
| 150-200 | 4.272 | 1.068 | 0.97650 | 13,699,644 |
| 200-250 | 2.008 | 1.142 | 0.90021 | 61,995,702 |
| 250-300 | 1.063 | 1.399 | 0.74095 | 221,149,984 |
| 300-350 | 0.552 | 1.976 | 0.40969 | 821,132,144 |
| 350-400 | 0.282 | 3.647 | 0.23986 | 3,137,860,490 |
| 400-450 | 0.152 | 6.599 | 0.07622 | 10,884,777,013 |
| 450-500 | 0.0795 | 12.490 | -0.05070 | 39,568,454,997 |

The Q=10 path projection assumes ideal Monte Carlo scaling
`Q(N)=Q(N0)*sqrt(N/N0)`. It is a lower-bound planning estimate once amplitude
and cosine also deteriorate. The current algorithm remains above aggregate
Q=10 through 150 a.u., but a one-billion-path budget cannot reach Q=10 over
the complete 500-a.u. interval without additional variance reduction.

## Files

- QM: `qm/ahm-qm-s10-n5.dat`
- Independent SepMB repeats: `run1/` and `run2/`
- Weighted mean: `combined/ahm-sepmb-s10-n5-2500000.dat`
- Window table: `analysis/long_window_metrics.csv`
- Q plot: `analysis/long_window_q.png`
- Full all-orbital plot: `analysis/v071_t500_all_orbitals_qm_vs_sepmb.png`
- Converged-region zoom: `analysis/v071_t500_zoom_0_200.png`
- Transition-region zoom: `analysis/v071_t500_zoom_0_300.png`
