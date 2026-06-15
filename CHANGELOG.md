# Changelog

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
