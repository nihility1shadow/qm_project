# 2500-au bounded-reference validation

## Scope

- System: 10 orbitals, 5 spinless electrons.
- Time grid: 0 to 2500 a.u., dt = 0.5 a.u., 5001 output rows.
- Bath: semi-elliptic hybridization, wc = 4 eV.
- Integrated coupling: eta = 0.04 eV^2 =
  5.4020507214803505e-5 Ha^2.
- Molecular energy: O 2p value, delE = -17.194874968839816 eV.
- Fermi energy: -4.5 eV.
- Vibrational basis: 384 levels.
- Acceptance metric:

  Q = RMS(QM occupation - initial occupation) /
      RMS(estimate - QM).

  Q is evaluated separately for every 100-a.u. window and every active
  orbital 0, 5, 6, 7, 8, and 9. The strict requirement is that the minimum of
  these scores is greater than 10. Orbitals 1-4 are still plotted and included
  in the absolute-error check, but their nearly zero physical signal makes a
  relative Q threshold ill-conditioned.

## Result

| Reference | Electronic states | Global active Q | Minimum active-window Q | Minimum active-orbital/window Q | Maximum orbital error | Reference propagation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| D2/F384 | 186 | 90.9941 | 30.2927 | 24.9959 | 2.37621e-6 | 7.08 s |
| D3/F384 | 246 | 90.9944 | 30.2928 | 24.9959 | 2.37619e-6 | 8.05 s |
| D4/F384 | 252 (full electronic space) | 90.9944 | 30.2928 | 24.9959 | 2.37619e-6 | 9.57 s |

All three references pass. For D3/F384, the worst case is orbital 8 in the
2400-2500 a.u. window, with Q = 24.9959174. D3 and full-electronic-space D4
differ by at most 6.66e-16 in orbital occupation; the remaining discrepancy
against grid QM is therefore not an electronic-distance truncation error.

The D3 job used 64 MPI ranks. Its maximum peak RSS per rank was 14,292 KiB.
The sum of individual rank peaks was 792,772 KiB (about 0.756 GiB), which is
not guaranteed to be simultaneous resident memory. The full scheduled job
took 54 s, including 42.53 s spent on a deliberately tiny 100-path residual
diagnostic. The bounded reference itself took 8.05 s. The one-rank grid-QM
comparison took 2,041 s (34 min 1 s).

## Final deployment regression

Server job 648262 tested the complete five-file source set after the cloud
matrix-vector compatibility fix. It used D3/F384 at 100 a.u. on 64 MPI
ranks. The bounded reference took 0.308 s; the complete computation took
0.421 s, and the maximum peak RSS per rank was 11,148 KiB.

All 201 by 14 values in `reference-observables.dat` are bit-for-bit identical
to the pre-deployment D3 regression from job 648246. The maximum absolute
difference and RMS observable difference are both zero. This confirms that
replacing the unavailable `CMatrix::vectorply` call with an explicit
matrix-vector loop fixed cloud compilation without changing the evolution.

## Algorithm interpretation

The reusable bounded propagator is analogous to the central idea of
Inchworm: use already resolved propagator information as a bold/control
quantity instead of resampling the full long-time history. This implementation
is not a complete temporal Inchworm expansion; it reuses a bounded
many-electron determinant subspace and samples only paths that leave it.

A bias-free n mod 4 phase-sector stratification was also tested as a
Modified-Filinov-inspired variance reduction. At 100 a.u. and 200,000 paths:

| Estimator | Global Q | Second-half Q | Sampling time |
| --- | ---: | ---: | ---: |
| Existing parity-CDF estimator | 9.5316 | 7.7871 | 320.53 s |
| Backward n mod 4 sectors | 3.6121 | 2.4987 | about 320 s |
| Forward and backward n mod 4 sectors | 3.0467 | 2.5458 | 328.7 s |

The existing shifted parity CDF already creates useful negative correlation
between opposite phase sectors. Splitting those sectors destroys that
correlation, so the experiment was rejected and the production source was
restored. A direct Filinov damping factor was not added because it would alter
the estimator rather than remain unbiased.

## Important limitation

The strict 2500-a.u. pass belongs to the deterministic D3/F384 bounded
reference. It is not evidence that the pure Poisson residual by itself reaches
Q > 10 at 2500 a.u. The 100-path residual in the probe is intentionally far
too small and is not used in the accepted result. The reference file itself
is the accepted single-run output.

For larger orbital/electron spaces, D3 may represent a much smaller fraction
of the Hilbert space. The same speed and accuracy cannot be claimed without a
new D2/D3/D4 convergence check and, where needed, a demonstrably beneficial
Poisson correction.

## Reproduction

Run the local analysis after the four job folders have been downloaded:

    python experiments/paper_semicircle_s10n5_t2500_20260911/analyze_reference_t2500.py

Outputs:

- validation.json: complete per-window and per-orbital metrics.
- summary.csv: compact acceptance table.
- reference_qm_all_orbitals.png: all ten occupation changes.
- reference_qm_quality.png: weakest active-orbital Q and maximum error.

## Rollback

- Local pre-experiment source:
  archive/phase-sector-before-20260911/ahm-mb-sep.cpp
- Rejected phase experiment:
  archive/phase-sector-before-20260911/ahm-mb-sep-failed-mod4.cpp
- Server pre-experiment source:
  /data/home/yd101802/yd101802/nonadia/_rollback_phase_sector_20260911/
- Server pre-deployment five-file source:
  /data/home/yd101802/yd101802/nonadia/_rollback_stage_2500_20260911/
- Restored server source SHA-256:
  86a1a407b2db495c4eaae8002a43dd54cbecacaea4d53556c20c490cc94ceb0c
