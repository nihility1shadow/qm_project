# One-run 100-million-trajectory result

## Run

- Job: `636249`
- Elapsed time: 12:26
- MPI ranks: 64
- Orbitals/electrons: `Norb=4`, `Nel=2`
- Trajectories: 100,000,000
- Parameters: `wc=0.25 eV`, `eta=6e-5`, `AHM_DELE_EV=-2.5`
- Time range: 0 to 100 a.u., 201 measured points at 0.5 a.u. intervals
- Maximum particle-number error: `1.64e-10`

## Interpretation

A single run cannot estimate repeat-to-repeat strict Q. The 100-million result was
therefore compared with the mean of two independent 30-million runs. For each
10 a.u. window and orbital, the cross-consistency metric is

`Q_cross = reference signal RMS / RMS(100m run - 30m mean)`.

The active orbitals are 0, 2, and 3. Their minimum cross-Q by window is:

| Window (a.u.) | Minimum active-orbital Q_cross |
|---:|---:|
| 0-10 | 823.13 |
| 10-20 | 254.69 |
| 20-30 | 83.08 |
| 30-40 | 79.34 |
| 40-50 | 39.35 |
| 50-60 | 27.45 |
| 60-70 | 18.69 |
| 70-80 | 15.87 |
| 80-90 | 13.22 |
| 90-100 | 11.11 |

Thus the single 100-million curve agrees with the two-run 30-million reference
through 100 a.u. at `Q_cross > 10`.

If Monte Carlo square-root scaling is extrapolated from the measured 30-million
strict result `Q=9.49`, a 100-million run is expected to have strict
`Q approximately 17.33`. This value is a prediction, not a measured strict Q,
because only one 100-million run was performed.

## Noise consistency

For independent Monte Carlo estimates, the expected RMS difference between a
100-million run and the mean of two 30-million runs, expressed in units of one
30-million-run standard deviation, is

`sqrt(1/2 + 30/100) = 0.894`.

The measured active-orbital mean ratio in the final 90-100 a.u. window is 0.996.
Given that the reference standard deviation is estimated from only two repeats,
this is consistent with ordinary sampling variation and does not indicate a
late-time systematic drift.

## Files

- Raw data: `cloud-runs/636249/ahm-sepmb-s4-n2-100000000.dat`
- Window metrics: `plots/single_vs_reference_windows.csv`
- All-orbital comparison: `plots/all_orbitals_30m_mean_vs_100m.png`
