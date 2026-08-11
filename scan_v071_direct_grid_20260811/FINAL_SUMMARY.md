# SepMB v0.71 direct-grid sampled-backward milestone

## Fixed physical configuration

- `Norb=10`, `Nel=5`
- `wc=3.367556 eV`
- `eta=3.2920651000724556e-5`
- `delE=-17.194874968839816 eV` (oxygen outermost-orbital setting)
- `dt=0.5 a.u.`, `nstep=400`, `tmax=200 a.u.`
- all 401 output times are evaluated directly (`SEP_MB_MEASURE_STRIDE=1`)

## Root cause of the performance problem

v0.64 propagated a dense backward vector over all `Nhs=C(10,5)=252`
many-electron determinants for every backward replica and every measurement
time.  The same full sparse-matrix/vector propagation was repeated many times,
even though one Monte Carlo projection only needs an unbiased estimate of the
backward amplitude at the determinant selected by the forward branch.

For an exact backward jump

```
v'(s') = sum_s J(s',s) v(s),
```

v0.71 samples one allowed target `s'` with probability `p(s'|s)` and carries

```
w' = w * J(s',s) / p(s'|s).
```

The target is uniform over the allowed determinant transitions.  If the degree
is `d`, the sampled multiplier is therefore `d*J(s',s)`.  Its expectation is
the exact matrix-vector result, so this changes variance and cost, not the
Hamiltonian or the expected observable.

The dominant zero- and one-jump backward sectors are still summed exactly with
a sparse active-state list (`SEP_MB_EXACT_BACK_JUMPS=1`).  This is a
Rao-Blackwell step: cheap low-order branches remain deterministic, while only
the expensive multi-jump orbital branching is sampled.

## Source changes

- `ahm-mb-sep.cpp:270-300`: add sampled-backward controls and production
  defaults (`B=256`, forward stratification on, exact orbital propagation on,
  single-jump-time stratification on, sampled backward orbitals on, exact
  backward jumps through order one).
- `ahm-mb-sep.cpp:431-485`: precompute determinant targets and fermionic sign /
  importance multipliers for every allowed center-to-bath or bath-to-center
  jump.
- `ahm-mb-sep.cpp:550-688`: use sparse exact propagation for low jump counts and
  a scalar determinant branch with inverse-probability weight for larger jump
  counts.
- `ahm-mb-sep.cpp:1040-1060`: write the active algorithm and every sampling
  switch into the data header.
- The backward-replica range is now `1..1024`; a fixed-budget screen selected
  `B=256` for this 10-orbital, 5-electron, 200-a.u. case.

## Fair 401-point benchmark

All values below use 64 MPI ranks and the same QM reference.

| method | forward paths | B | runtime | global Q | RMS error | cosine |
|---|---:|---:|---:|---:|---:|---:|
| QM reference | - | - | 166 s | - | - | - |
| v0.64 dense backward | 200,000 | 64 | 527 s | 1.1538 | 1.4594e-5 | 0.7618 |
| v0.71 sampled backward | 100,000 | 256 | 138 s | 1.5781 | 1.0670e-5 | 0.8551 |

v0.71 is 3.82 times faster than the retained v0.64 benchmark and 1.20 times
faster than this small 252-state QM calculation, while its RMS error is 26.9%
lower than v0.64 in this run.

Window quality for the final v0.71 run:

| time window (a.u.) | Q |
|---|---:|
| 0-50 | 12.48 |
| 50-100 | 5.437 |
| 100-150 | 2.090 |
| 150-200 | 0.8818 |

The late-time sign/phase problem is reduced but not removed.  This example is
well converged at early times and improves the middle windows, but it is not a
claim of `Q>10` through 200 a.u.

## Default server command

No new algorithm environment variables are required:

```bash
sbatch --ntasks=64 \
  --export=ALL,AHM_WC_EV=3.367556,AHM_ETA=3.2920651000724556e-5,AHM_DELE_EV=-17.194874968839816,AHM_NSTEP=400 \
  run_mpi_cloud.slurm 100000 10 5 1
```

## Server rollback

Restore the original direct-grid v0.64 implementation:

```bash
cp ahm-mb-sep-v064.cpp ahm-mb-sep.cpp
cp na_mpi_v064.out na_mpi_cloud.out
```

Restore the faster interpolated-grid v0.70 checkpoint:

```bash
cp ahm-mb-sep-v070.cpp ahm-mb-sep.cpp
cp na_mpi_v070.out na_mpi_cloud.out
```

Re-activate this milestone:

```bash
cp ahm-mb-sep-v071.cpp ahm-mb-sep.cpp
cp na_mpi_v071.out na_mpi_cloud.out
```

Server SHA-256 identifiers:

- v0.71 source: `63179a2b501f2bf2de1b2dc84219538825aeba49b52fa1a9d4a6904685483e67`
- v0.71 executable: `c63ac360b5cb06d08241bfdc464e19323a308f1c2e7c4caf68829dc9fd87ed38`
- v0.64 source: `ae3b42e26c19454a118ff12a562ee07d6c8c25ec1dda86047e6b289641d0d43a`
- v0.64 executable: `61aff826d52b116a6933d2692f1d5b4758bf6251ad405a208354404210478f5c`
