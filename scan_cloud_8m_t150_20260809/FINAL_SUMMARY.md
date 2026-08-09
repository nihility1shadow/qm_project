# 8M-trajectory, 150 a.u. parameter search

## Fixed conditions

- `Norb = 10`, `Nel = 5`
- `ntraj = 8,000,000` per independent run
- `t = 0..150 a.u.` with real measurements every `0.5 a.u.`
- `AHM_DELE_EV = -2.5 eV`
- strict score: minimum active-orbital Q over every complete 10 a.u. window through 130 a.u.
- pass requirement: strict Q greater than 4
- each search round retained the top 30%, rounded upward

## Search history

| Stage | Candidate pool | Repeats | Best strict Q | Best parameters |
|---|---:|---:|---:|---|
| Round 1 | 12 | 3 | 2.105 | `wc=0.20 eV`, `eta=5.0e-5` |
| Round 2 | 16 | 3 | 2.151 | `wc=0.20 eV`, `eta=4.0e-5` |
| Round 3 | 17 | 3 | 2.230 | `wc=0.20 eV`, `eta=3.5e-5` |
| Round 4 | 18 | 3 | 2.335 | `wc=0.12 eV`, `eta=2.0e-5` |
| Final validation | 6 | 6 | 2.323 | `wc=0.12 eV`, `eta=2.8e-5` |

## Final retained top 30%

1. `wc=0.12 eV`, `eta=2.8e-5`, strict Q = `2.3233`
2. `wc=0.20 eV`, `eta=3.75e-5`, strict Q = `2.2619`

Neither candidate passes Q > 4 through 130 a.u. The validated winner has Q =
4.2869 in 80-90 a.u., 3.4631 in 90-100 a.u., and 2.3233 in 120-130 a.u.

Using the observed Monte Carlo scaling `Q proportional to sqrt(ntraj)`, the direct
estimate for the winner is

`ntraj(Q=4) = 8,000,000 * (4 / 2.3233)^2 = 23.7 million`.

A 32-million-trajectory validation is the conservative next test. This remains
well below the stated one-billion-trajectory ceiling.

## Physical caveat

The numerical winner is an ultranarrow band. With `EF=-4.5 eV` and `wc=0.12 eV`,
the discretized bath lies approximately in `[-4.56, -4.44] eV`, so all bath
orbitals are close to the Fermi level. It maximizes the present convergence score,
but it should not be interpreted as evidence that a broad physical reservoir has
only two active levels.
