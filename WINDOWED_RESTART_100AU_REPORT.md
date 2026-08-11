# 4 x 25 a.u. window-restart experiment

## Objective

Test whether the accurately converged short-time Poisson state can be used as
the initial state of the next time segment.  The test covers 100 a.u. as four
25 a.u. windows and compares against both exact QM and the active continuous
SepMB v0.71 implementation.

The quality metric is

\[
Q = \frac{\operatorname{RMS}\left[n_j^{\rm QM}(t)-n_j(0)\right]}
         {\operatorname{RMS}\left[n_j^{\rm method}(t)-n_j^{\rm QM}(t)\right]},
\]

where the RMS includes all ten orbitals and all time points in a 25 a.u.
window.  Larger values are better, and the target is Q > 10.

## Fixed parameters

- `Norb = 10`, `Nel = 5`
- `dt = 0.5 a.u.`, `tmax = 100 a.u.`
- `wc = 3.389460 eV`
- `eta = 2.961e-5`
- `delE = -17.194874968839816 eV`
- `ntraj = 100000`
- Continuous v0.71: 64 conditional backward replicas
- Window v0.75: two independent 100000-particle ket/bra ensembles, eight
  class-conditioned pairings, four 50-step windows
- Final window result does not merge different coherent-state alpha values.

## Implementation tested

At a window boundary, every particle carries the many-electron amplitude
vector, coherent-state parameter alpha, complex path weight, and excitation
class.  Systematic importance resampling creates the ensemble used as the
initial state of the next window.  Orbital occupations alone are never used as
the restart state.

Several corrections were needed during the experiment:

1. Randomly permute or condition ket-bra pairings instead of pairing identical
   stratification indices.
2. Balance the jump-count low-discrepancy sequence within every MPI rank.
3. Preserve complex phases during boundary resampling.
4. Reject alpha-bin merging at widths 0.05 and 0.10 because its approximation
   error was larger than its variance reduction.

All window code is disabled by default.  It is enabled only when
`SEP_MB_WINDOW_STEPS` is positive, so normal execution retains v0.71 behavior.

## Accuracy

| Time window (a.u.) | Continuous v0.71 Q | 4 x 25 v0.75 Q |
|---:|---:|---:|
| 0-25 | 34.623 | 9.493 |
| 25-50 | 9.290 | 0.774 |
| 50-75 | 3.498 | 0.317 |
| 75-100 | 1.702 | 0.204 |

The four-window method conserves the total electron count to approximately
`3.2e-12`, but conservation alone does not imply accurate orbital dynamics.
Its noise grows after every restart boundary and clearly exceeds v0.71 after
25 a.u.

Assuming ideal independent Monte Carlo scaling, the last window would require
approximately 240 million forward trajectories to reach Q = 10.  This is an
optimistic estimate because boundary errors are inherited and correlated, so
the observed error is not purely independent `1/sqrt(N)` noise.

## Timing

| Method | MPI ranks | Wall time (s) | User CPU time (s) |
|---|---:|---:|---:|
| Exact QM | 1 | 83.05 | 81.68 |
| Continuous v0.71 | 128 | 8.25 | 377.76 |
| 4 x 25 v0.75 | 128 | 5.43 | 193.15 |

The window prototype is about 34% faster in wall time than v0.71 for this
100000-trajectory test, but the accuracy loss is much larger than the speed
gain.

## Conclusion

The physical identity

\[
|\Psi(t+\Delta t)\rangle = U(t+\Delta t,t)|\Psi(t)\rangle
\]

is valid, so exact state segmentation is valid.  The tested finite-particle
restart is not equivalent to exact state segmentation: statistical phase
cancellation that has not converged at one boundary becomes part of the next
window's initial state.  Resampling then propagates and amplifies that boundary
error.  The current v0.71 estimator remains the recommended active version.

A useful future window method needs to compose a converged short-time
propagator or influence kernel (an inchworm-style construction), rather than
restart from a finite resampled particle cloud.  Such a method must retain the
bra-ket correlations and electron-vibration memory across boundaries.

## Rollback

- Server active source and executable remain `ahm-mb-sep.cpp` and
  `na_mpi_cloud.out`, both identical to their v0.71 backups.
- Server rollback files are `ahm-mb-sep-v071.cpp` and `na_mpi_v071.out`.
- Window v0.75 is an isolated candidate (`ahm-mb-sep-v075-candidate.cpp` and
  `na_mpi_v075_test.out`) and is not deployed as the active program.
