# Parameter optimization result

## Constraints

- Bath cutoff: `wc <= 10 eV`.
- Required interval: every 10 a.u. segment through `60 a.u.` must have `Q >= 2`.
- Particle-number maximum error must remain below `1e-8`.
- Final validation uses 8,000,000 stochastic trajectories per run, 64 MPI ranks,
  four independent repeats, and output times `0, 0.5, ..., 75`.

## Selected parameter

- `wc = 4.5 eV`
- `eta = 1.20e-4`
- Minimum `Q` through 60 a.u.: `2.5892514198`
- `Q(40-50) = 2.5892514198`
- `Q(50-60) = 2.9231509638`
- `Q(60-70) = 1.5947147488`
- `Q(70-75) = 1.1151381353`
- Particle-number maximum error: `7.9602102687e-11`

The required interval is therefore valid through 60 a.u.  The data after 60
a.u. are retained for trend inspection but are not classified as converged.

## Runner-up

- `wc = 4.4 eV`
- `eta = 1.15e-4`
- Minimum `Q` through 60 a.u.: `2.3892598484`
- `Q(50-60) = 2.4427577402`

## Screening conclusion

The high-cutoff branch (`wc` around 8.2 eV) ranked first at 50 a.u. but fell
to `Q(50-60) = 1.5784603869` at 75 a.u.  The lower-cutoff, higher-eta branch
is more stable at longer times.  Interpolation between 4.2 and 4.6 eV located
the current optimum at 4.5 eV and `eta = 1.20e-4`.

Numeric convergence and plotting are separate.  `rank_sci_q.py` writes only
CSV/JSON metrics; `plot_scan_results.py` creates PNG plots on demand.  The C++
flat-tail analyzer no longer writes HTML unless `--html-report` is requested.
