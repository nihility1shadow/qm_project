# v1.12 compact bounded reference: results and remaining accuracy gap

## Scope

Preserves the physical Hamiltonian, fermionic signs, Kondo/star Poisson
paths and the deterministic-reference plus sampled-residual estimator.
The user accepts 1-2 days per run, requires GitHub rollback stages and plots.
Context and physical parameters are in V111_DISTRIBUTED_REFERENCE.md.

For path-local systems with at most 64 orbitals, reference determinants now
use 64-bit masks and hash lookup instead of permanent set/map trees.
The greater-than-64-orbital set implementation remains available.
Reference traversal order is unchanged. This is bounded reference storage,
not an allocation of the full many-electron basis.

A backward path that cannot leave the reference has an exactly zero residual
when the forward path also stays in it. Such a path skips coherent propagation
only after all random draws, preserving the sampling sequence. Set
SEP_MB_SKIP_ZERO_REFERENCE=0 to disable this work elimination.

## Correctness checks

Run `python tests/test_reference_masks.py`: 652 graph/order cases pass,
including non-half filling, an initially empty impurity and orbital bit 63.
Remote MPI build uses GCC 12.1.0 and OpenMPI 4.1.8.

All following observable-file comparisons have maximum absolute difference 0:

| Jobs | Comparison |
| --- | --- |
| 647607 / 647618 | v1.11 / v1.12, 10/5, 100 a.u., D3 |
| 647608 / 647619 | v1.11 / v1.12, 30/15, 10 a.u., D2 |
| 647620 / 647621 | zero-residual skip off/on, 10/5, 500 a.u., D1 |
| 647634 / 647635 | continuous-H reference MPI off/on, 8/3, 10 a.u. |
| 647636 / 647637 | H0/H1 reference MPI off/on, 8/3, 10 a.u. |

The existing uneven 386-level/7-rank test differs by only 2.22e-15.
The 10/5 comparison to same-parameter grid QM is retained from v1.11:
active Q=2.50e6 over 100 a.u., but its reference is 246 of 252 states.
This small-system result does not establish 30/15 accuracy.

## Resource measurements

30 orbitals / 15 electrons, D2=52,656 states, 384 oscillator levels, 64 ranks:

| Test | Wall time | Maximum rank peak | Sum of rank peaks |
| --- | ---: | ---: | ---: |
| v1.11, 10 a.u., 100 paths | 34.17 s | 144.51 MiB | 8.92 GiB |
| v1.12, 10 a.u., 100 paths | 26.10 s | 41.07 MiB | 2.46 GiB |
| v1.12, 2500 a.u., 100,000 paths | 2h 05m 15s | 44.19 MiB | 2.65 GiB |

The long v1.12 run is job 647622. Two same-settings v1.11 seeds (647613,
647614) took 2h26m and 2h28m. Since the seeds differ, the short matched-seed
pair is the controlled comparison. Individual rank peaks are not simultaneous
RSS; their sum is a conservative resource summary, not a simultaneous sample.

A D3 capacity test (647624) used 715,136 reference states, 10 a.u., 100 paths:
369.30 s, 444.80 MiB/rank maximum, 27.68 GiB sum on 64 ranks. A linear
reference-time projection to 2500 a.u. is approximately 25.5 hours, plus
sampling. This is a projection, not a completed long D3 validation.
30/15 full determinant count is 155,117,520; D3 covers 0.461% of that space.

## Completed long-time validation: D2 does not yet meet the accuracy goal

Three independent seeds (647613, 647614, 647622), each 100,000 paths,
2500 a.u., 52,656 reference states. Version equivalence is established above.
Q_repeat = RMS(mean occupation change) / RMS(single-run standard deviation),
computed with ddof=1 across the three runs in every 100-a.u. window.
This is NOT Q against exact QM. It cannot detect a shared deterministic bias,
and three runs alone cannot establish the tail behavior of rare path weights.
Active orbitals are 0 and 15-29; all 30 are plotted and particle number checked.

Over the last window [2400,2500], combined active Q_repeat=4.446 and the
weakest active orbital Q_repeat=0.676. Every active orbital exceeds 10 only
through the [1900,2000] window. Thus even the diagnostic Q=10 threshold fails
at late times. Do not advertise this version as high-fit over all 2500 a.u.
Maximum particle-number error is 5.16e-14. Some estimates become slightly
negative (minimum -5.87e-7); no clipping or smoothing is applied.

The matched-seed 384/512 oscillator test (647613/647623) differs by at most
2.816e-11 per orbital, much less than late sampling noise. It checks oscillator
truncation at D2, not reference-space convergence or absolute QM accuracy.
Increasing paths alone to raise weakest Q_repeat to 10 is roughly a 219-fold
cost projection, assuming ideal inverse-square-root noise scaling. Enlarging
the bounded reference is the next candidate within the accepted time budget.

## Reproduction and rollback

`python analyze_v111.py` and `python analyze_v112_long.py` regenerate checked
JSON reports and static PNGs under scan_v111_distributed_20260905/figures.
Selected raw observables, configuration, logs, binary hashes and timing are
tracked. Sampler diagnostic files are excluded from fitting comparisons.
`run_poisson_lowmem.slurm` selects the separate v1.12 candidate binary; the
server's historical default executable is not replaced. Its arguments are:

```
steps paths orbitals electrons reference_distance max_reference_states mpi seed
```

Example: 64 ranks, 384 levels, D2, 2500 a.u.:

```bash
sbatch --ntasks=64 --time=24:00:00 run_poisson_lowmem.slurm 5000 100000 30 15 2 65536 1 111103
```

Rollback stages: v1.10-path-local-large-space, v1.11-mpi-reference,
v1.12-compact-reference. Switch to a tag in an isolated checkout to preserve
unrelated local work. High-order long validation is the next stage.
