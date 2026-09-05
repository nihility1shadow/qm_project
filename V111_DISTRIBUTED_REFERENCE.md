# v1.11: distributed bounded reference for Poisson dynamics

## Inherited context

Read the accessible project task history: 269 turns in the main debugging task,
6 in the short-time-kernel discussion, 2 in the oxygen-orbital discussion and
4 in the archived convergence-analysis task.

Current constraints: preserve the Anderson–Holstein Hamiltonian and Kondo/star
Poisson paths; do not allocate the complete determinant basis; judge fitting on
orbital changes rather than the near-unit occupation baseline. Keep oxygen
delE=-0.6319 Hartree and wc <= 10 eV. Show every orbital, actual wall time and
peak memory. The current accepted runtime budget is one to two days per run.
Commit tested stages to GitHub and preserve rollback.

The previous 30/15 run used 1696 reference states, 384 oscillator levels and
100,000 paths over 2500 a.u.; 51m58s on 64 ranks, maximum rank RSS 46.4 MiB.
Its late-time noise was unresolved. No full-QM truth exists for that run.
The small 10/5 success does not establish large-system accuracy.

## Change and limits

Opt in with SEP_MB_REFERENCE_MPI=1. Each MPI rank owns consecutive oscillator
levels and three local complex work vectors. Only two oscillator rows at each
boundary are exchanged; electronic hopping remains local. Global oscillator
indices are used for coefficients and truncation. Observables are reduced to
the master rank. The Poisson sampler and its random-number consumption are
unchanged. SEP_MB_REFERENCE_MPI=0 retains the serial reference path.

This is a deterministic-reference plus Poisson-residual method, not pure
Poisson sampling. Distributing the reference changes its cost, not the
estimator. Enlarging the bounded reference can reduce the residual variance;
the long runs must establish how much. It does not solve the general
real-time sign problem or guarantee accuracy for arbitrary Hamiltonians.

Reference storage per rank is approximately
16 * Nref * (3 * local_Nfock + 4) bytes for MPI wavefunction/work/halo arrays,
plus replicated determinant and sparse-transition metadata. The process peak
therefore does not scale solely with the wavefunction formula. At least two
oscillator levels per rank are required. MPI remains opt-in.

All individual rank peaks are summed in the log. This sum is an upper bound on
simultaneous total RSS, not a measurement of simultaneous memory consumption.
External time.txt measures the whole launched simulation.

## Validation

- Same seed, 4 ranks, 10/5, 100 a.u., 1000 paths: old v1.10 versus v1.11 MPI off,
  and MPI off versus on: maximum output difference 0.
- 386 oscillator levels on 7 ranks (unequal partitions), 10/5, 10 a.u.:
  maximum output difference 2.22e-15 (roundoff).
- 30/15, 10 a.u., 100 paths, 64 ranks: 52656 reference states at distance 2,
  147980 KiB maximum rank peak RSS; sum of rank peaks 9353508 KiB.
- See stage1_validation.json and figures for measured timing and same-parameter
  small-system grid-QM comparison.

The initial old-binary run failed to locate GSL before simulation began.
run_v111.slurm now supplies the project dependency directory explicitly;
the rerun passed. That failed launch is not included in performance claims.

Long pilot jobs 647613 and 647614: 30/15, 2500 a.u., 100000 paths each,
64 ranks each, independent seeds 111101/111102, distance 2 / 52656 reference
states / 384 oscillator levels. Results are pending at this stage.

## Reproduce and rollback

run_v111.slurm arguments:
nstep ntraj Norb Nel reference_distance max_reference_states reference_mpi seed

Example: sbatch --ntasks=64 --time=24:00:00 run_v111.slurm 5000 100000 30 15 2 65536 1 111101

Run analyze_v111.py locally to regenerate the stage figures.
The default server executable is not replaced. Candidate executable:
na_mpi_v111_distributed.out. Previous executable:
na_mpi_v110_path_local.out. Previous Git checkpoint: 7c44668 and
v1.10-path-local-large-space. A separate checkout of either tag is the safest
rollback when the working tree contains unrelated changes.

## Method context

The general deterministic-core/stochastic-remainder idea is related to
[Petruzielo et al., Semistochastic Projector Monte Carlo Method](https://arxiv.org/abs/1207.6138)
and [Blunt et al., semi-stochastic FCIQMC developments](https://arxiv.org/abs/1502.04847).
Those papers concern projector methods; they do not demonstrate the accuracy
or time scaling of this project's real-time Poisson estimator.
