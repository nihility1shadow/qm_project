# Changelog

## v0.11.0 - Match SepMB output time grid to QM data

Function:
- Make separated many-body output use the full simulation time grid so it can be compared row-by-row with `ahm-qm-s10-n5.dat`.

Changes:
- Set `nwf = nstep` and `nbloc = 1` in `SepMBpoisson` instead of compressing output to 200 blocks.
- Add `ahm-mb-sep-v0.31.cpp` as the full-time-grid snapshot.
- Update the `#PATCH_CHECK` marker to `v0.31 full-time-grid.

## v0.10.0 - Prepare server overwrite source

Function:
- Provide the corrected separated many-body source under the original server filename while preserving the pre-v0.30 file for rollback.

Changes:
- Add `ahm-mb-sep.cpp` as a copy of the verified `ahm-mb-sep-v0.30.cpp` so it can directly overwrite the server file.
- Add `archive/ahm-mb-sep-before-v030.cpp` containing the previous desktop/source version for explicit rollback reference.
- Keep `ahm-mb-sep-v0.30.cpp` as a named fixed snapshot.

## v0.9.0 - Fix separated many-body normalization

Function:
- Add a corrected `SepMBpoisson` source version and verified output for the separated many-body estimator.

Changes:
- Add `ahm-mb-sep-v0.30.cpp` as a separate fixed version without overwriting the original file.
- Keep Poisson compensation in `sclf` without embedding `1/ntraj`.
- Apply the Monte Carlo average exactly once in the forward weight passed to `csproj`.
- Keep the backward final-interval compensation as `sclf[j-offset]`; using full `sclf[j]` overcompensates and grows the trace.
- Normalize printed observables by `trace = prb[t][0] / Nel`, then print the conserved electron count as `Nel`.
- Fix `jc` and MPI reductions to use `nstep + 1` entries.
- Set `Jmax` to `nstep + 1` so sampled jump arrays cannot overflow.
- Fix the `nj_min` parity expression with explicit parentheses.
- Add a `#PATCH_CHECK` marker to generated `ahm-sepmb` files so stale binaries are obvious.
- Add verified output `ahm-sepmb-s10-n5-50000-v030.dat`.

## v0.8.0 - Harden array allocation helpers

Function:
- Make multi-dimensional allocation helpers fail predictably when memory allocation fails.

Changes:
- Check `array2d` row-pointer allocation before dereferencing it.
- Report separate allocation failures for pointer layers and data blocks.
- Add allocation checks and cleanup for `array3d` and `array4d`.

## v0.7.0 - Repair Kondo DP table construction

Function:
- Make the Kondo path sampler DP probabilities safe and consistent with the documented recurrence.

Changes:
- Remove incorrect guards that skipped valid transitions from adjacent distance states.
- Fix `Qtd` boundary checks to avoid negative and out-of-range distance indexes.
- Handle degenerate `M == 0` or `N-M == 0` cases before division.
- Initialize sampler table pointers before building DP tables.
- Add an include guard for `Kondo-path-sampler.h`.
- Disable the malformed experimental `compute_dp_table_sum` body and replace it with a safe wrapper.

## v0.6.0 - Guard FFT plan destruction

Function:
- Make default-constructed `QMFFT1D` objects safe to destroy.

Changes:
- Initialize `fftpln_for` and `fftpln_bck` to `NULL` in the default constructor.
- Destroy FFTW plans only when the plan handles are non-null.

## v0.5.0 - Restore Johnson path sampler

Function:
- Provide a working Johnson graph path sampler used by the Kondo path sampler.

Changes:
- Replace the placeholder `Johnson-path-sampler.cpp` class declaration copy with function implementations.
- Add the Johnson distance DP table and conditional path sampling logic.
- Implement `get_random_element`, `print_path`, and `verify_path`.
- Add a header guard and release the allocated `Ptd` table in the destructor.
- Define the shared `log_table` storage used by the sampler.

## v0.4.0 - Fix 1D reallocation helper

Function:
- Restore `realloc1d` so code using dynamic 1D reallocation can compile and allocate correctly.

Changes:
- Replace the invalid `malloc(ptr, size)` call with `malloc(size)` in `yyy_inlines.h`.

## v0.3.0 - Version notes

Function:
- Add project version notes so each pushed version records its purpose and changes.

Changes:
- Document existing versions and the Monte Carlo averaging bug fix.

## v0.2.0 - Accumulate Monte Carlo wavefunction samples

Function:
- Fix the multi-trajectory Monte Carlo wavefunction average in `AHM::MBpoisson`.

Changes:
- Generate each trajectory contribution into a temporary buffer before adding it to `wf`.
- Keep the existing `1/ntraj` scaling while ensuring all trajectory contributions are accumulated.
- Allocate `jc` with `nstep + 1` entries and reduce/write the same length.
- Free 3D wavefunction arrays with `free3d`.
- Release auxiliary arrays before returning.

## v0.1.0 - Initial project import

Function:
- Import the quantum mechanics simulation source files and reference/output data.

Changes:
- Add AHM, QHO, FFT, quantum jump, Johnson path sampler, and Kondo path sampler source/header files.
- Add expected data file `ahm-qm-s10-n5.dat`.
- Add simulated data file `text 1000 10 5 (0).dat`.
