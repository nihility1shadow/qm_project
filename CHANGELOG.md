# Changelog

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
