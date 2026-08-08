# Strict per-orbital parameter validation

## Model and scan constraints

- Model: multi-electron Anderson-Holstein/Newns-Anderson calculation.
- `Norb = 10`, `Nel = 5`; these are the finite model-space size and total
  electron count, not an instruction to select a five-electron atom.
- The current molecular/local level is unchanged at the server default:
  `delE = -17.194874968839816 eV` (O atomic 2p HF-Roothaan value).
- Bath cutoff constraint: `0 < wc <= 10 eV`.
- Final validation: `32,000,000` trajectories per run, `128` MPI ranks,
  four independent repeats, `dt = 0.5 a.u.`, and output through `75 a.u.`.

## Strict convergence rule

The old pooled signal/noise score was discarded because a large-amplitude
orbital could hide a failed weak orbital. Active orbitals are classified from
their RMS occupation change over `0-40 a.u.` using a 5% relative threshold.
For this calculation they are orbitals `0, 5, 6, 7, 8, 9`.

For every active orbital and every 10 a.u. segment,

`Q_j = RMS_t[mean_r(delta n_j)] / RMS_t[SD_r(delta n_j)]`.

An interval passes only when every active orbital has `Q_j >= 2`, inactive
orbital noise stays below 5% of the strongest reference signal, all four data
files are hash-distinct, and particle number is conserved within `1e-8`.

## Selected parameter

- `wc = 5.0 eV`
- `eta = 2.40e-4`
- Active orbitals: `0, 5, 6, 7, 8, 9`
- Maximum particle-number error: `1.10685e-10`
- Maximum inactive-orbital repeat noise through 60 a.u.: `9.393e-6`
- Inactive-noise limit: `2.548e-5`

| Time segment (a.u.) | Minimum active-orbital Q |
|---|---:|
| 0-10 | 84.157 |
| 10-20 | 18.236 |
| 20-30 | 7.126 |
| 30-40 | 4.888 |
| 40-50 | 3.165 |
| 50-60 | 2.501 |
| 60-70 | 1.019 |
| 70-75 | 1.042 |

The result is therefore validated through `60 a.u.`. Data after 60 a.u. are
retained for visual trend inspection but are not classified as converged.

Each final run used 128 ranks and completed in 91-98 seconds. The four Slurm
jobs were `635121`, `635122`, `635123`, and `635124`.

## Screening conclusion

The previous choice `wc = 4.5 eV`, `eta = 1.20e-4` is valid only through
40 a.u. under the strict rule. A literature-motivated per-level coupling of
`0.01 a.u.` corresponds to `eta = 9e-4` for nine bath orbitals, but the
multi-electron implementation showed that repeat variance grew faster than
the useful signal at `eta = 6e-4` to `9e-4`. The useful numerical region was
instead narrowed to `wc = 4-5 eV`, `eta = 2.4e-4` to `3.6e-4`.

The physical parameter scale was informed by the O2/silver Newns-Anderson
example in the project manuscript and by the continuous-metal-bath model in
[Huang, Xu, and Zhou (JCP 474, 111771)](https://arxiv.org/abs/2206.02173).
The literature supplies the model structure and order of magnitude, while the
final parameter is selected by independent-repeat convergence in this
multi-electron code rather than copied from a one-electron calculation.

## Random-seed correction and rollback

Simultaneous Slurm jobs previously used only the startup second and MPI rank,
so independent repeats launched in the same second could be byte-identical.
The cloud entry now mixes microseconds, `SLURM_JOB_ID`, and MPI rank. Setting
`AHM_SEED` explicitly reproduces the same output exactly. The test used two
simultaneous jobs: automatic seeds produced different SHA-256 hashes, while
`AHM_SEED=123456789` produced identical hashes.

The pre-fix server entry and binary are stored at:

`/data/home/yd101802/yd101802/nonadia/_rollback_seed_20260809_0115`

The server retains only the four final data directories `635121-635124`.

## Plots

- `plots/strict_validation_q60_best_orbitals.png`
- `plots/strict_validation_q60_ranking.png`
