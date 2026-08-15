# SepMB v0.81 500 a.u. validation

## Result

The continuous (non-window-restarted) SepMB algorithm reaches the requested
quality target over 0-400 a.u. with 45,000,000 forward trajectories:

- minimum aggregate active-orbital Q: 22.733
- minimum individual active-orbital Q: 17.173
- all 50 a.u. windows through 500 a.u. remain above Q=10
- full-range amplitude ratio: 1.00419
- full-range cosine against matching QM: 0.999672
- maximum particle-number error: 1.918e-8

Active orbitals are 0 and 5-9. Orbitals 1-4 have nearly zero QM signal; their
absolute SepMB fluctuations are around 1e-8 and are therefore assessed by
absolute error rather than an unstable relative Q.

## Final configuration

| parameter | value |
|---|---:|
| Norb / Nel | 10 / 5 |
| delE | -17.194874968839816 eV |
| wc | 2.7 eV |
| eta | 2.0e-6 |
| dt / nstep / tmax | 0.5 / 1000 / 500 a.u. |
| forward trajectories | 45,000,000 (3 x 15,000,000) |
| backward replicas | 256 |
| exact backward jumps | 4 |
| real measurement points | 180 adaptive points |
| output points | 1001 (0, 0.5, ..., 500 a.u.) |

The adaptive schedule directly samples every 0.5 a.u. through 63.5 a.u. and
then reduces the frequency according to the smallest electronic gap. Applying
the same schedule to noise-free QM gives a minimum interpolation-only Q of
250.1 over 0-400 a.u.; Monte Carlo variance, not interpolation, remains the
dominant error.

## Window Q

| window (a.u.) | aggregate Q | minimum active-orbital Q |
|---|---:|---:|
| 0-50 | 338.670 | 329.230 |
| 50-100 | 154.500 | 79.067 |
| 100-150 | 102.990 | 79.167 |
| 150-200 | 67.532 | 55.597 |
| 200-250 | 48.210 | 37.397 |
| 250-300 | 79.335 | 42.391 |
| 300-350 | 27.061 | 18.300 |
| 350-400 | 22.733 | 17.173 |
| 400-450 | 33.256 | 22.137 |
| 450-500 | 24.709 | 18.997 |

## Algorithm decisions

- Kept the original continuous trajectory evolution; segmented state restarts
  were rejected because they do not preserve the path-history estimator.
- Kept v0.76 sparse jump-time generation and removed redundant buffer clearing.
- Set exact_back_jumps=4 and adaptive measurement as v0.81 defaults.
- Rejected exact-zero backward sampling: Q was unchanged while runtime doubled.
- Rejected hashed backward-time starts: fixed-seed Q decreased.
- Kept B=256: B=1024 improved per-path precision but had worse time-to-target.

## Runtime and rollback

The three 15,000,000-trajectory jobs used 128 MPI ranks on two nodes each and
finished concurrently in 39-40 minutes. The server active files are v0.81;
v0.71 source, executable, and run script remain available for rollback.
