# SepMB v0.64 t=200 sampling requirement

## Fixed physical and numerical setup

- `Norb = 10`, `Nel = 5`
- `delE = -17.194874968839816 eV` (oxygen outermost-orbital energy)
- `wc = 3.367556 eV` (satisfies `wc <= 10 eV`)
- `eta = 3.2920651000724556e-5`
- `dt = 0.5 a.u.`, `nstep = 400`, `tmax = 200 a.u.`
- Every output point is directly evaluated: `SEP_MB_MEASURE_STRIDE=1`
- Active algorithm: SepMB v0.64, exact orbital sum and 2D one-jump-time stratification
- Long-time budget choice: `SEP_MB_BACK_REPLICAS=64`

The Q requirement is evaluated on the dynamically active orbitals
`{0, 5, 6, 7, 8, 9}`. Orbitals 1--4 have nearly zero QM signal, so their
quality is judged with absolute leakage rather than a relative Q ratio.

## Measured validation

The final validation used `Ntraj=200000`, `B=64`, and 64 MPI tasks. Slurm job
`644990` completed successfully in 527 seconds. The matching one-process QM
job `644979` completed in 166 seconds.

| Time window (a.u.) | Measured Q | Ntraj estimated for Q=10 |
|---:|---:|---:|
| 0--50 | 14.4367 | 95,961 |
| 50--100 | 2.9491 | 2,299,629 |
| 100--150 | 1.6638 | 7,224,807 |
| 150--200 | 0.6173 | 52,490,966 |

The estimate uses the verified Monte Carlo scaling

```text
Q(N) = Q(N0) * sqrt(N/N0)
N_required = N0 * (Q_target/Q(N0))^2.
```

An independent comparison between 25,000 and 200,000 trajectories agrees
with the expected `sqrt(8)` improvement in all four windows to useful
accuracy. Therefore the scaling estimate is suitable for run planning.

## Recommended production sizes

For the existing requirement that the simulation spans 200 a.u. and remains
at `Q > 10` through 150 a.u.:

- mathematical lower estimate: about 7.23 million forward trajectories;
- recommended size with margin: 10 million forward trajectories, `B=64`;
- predicted window Q at 10 million: 102.1, 20.85, 11.77, and 4.36.

Thus 10 million trajectories should meet Q>10 through 150 a.u., but not in
the final 150--200 a.u. window. It remains below the 24-million trajectory
limit specified for the earlier t=200 requirement.

If Q>10 is required throughout the complete 0--200 a.u. interval:

- mathematical lower estimate: about 52.5 million forward trajectories;
- recommended size with margin: 60 million forward trajectories;
- predicted final-window Q at 60 million: about 10.69.

Here `Ntraj` is the forward-trajectory argument passed to the executable.
With `B=64`, the internal backward-projection budget is approximately
`64*Ntraj`.

## Runtime estimate

Using the measured 527 seconds for 200,000 trajectories on 64 cores:

| Target | One 64-core job | Four independent 64-core jobs in parallel |
|---|---:|---:|
| 10 million total | about 7 h 19 min | about 1 h 50 min (4 x 2.5 million) |
| 60 million total | about 43 h 55 min | about 10 h 59 min (4 x 15 million) |

Queue waiting time is not included. Independent jobs can be averaged locally
and are statistically equivalent to using the summed trajectory count, while
also providing repeat-to-repeat error bars.

## Local outputs

- `validation/n200k-b64-644990/ahm-sepmb-s10-n5-200000.dat`: final Poisson validation
- `qm/ahm-qm-s10-n5.dat`: same-parameter t=200 QM reference
- `validation-analysis/v064_t200_n200k_b64_vs_qm.png`: all-orbital comparison
- `validation-analysis/t200_sampling_requirement.csv`: machine-readable estimates
- `replica-analysis/`: B32/B64 repeat-allocation comparison
