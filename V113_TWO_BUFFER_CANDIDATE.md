# v1.13: two-buffer reference implementation

The original candidate is preserved at tag v1.13-csr-candidate. Full cloud
regressions and the small-system 4000-a.u. comparison have now completed.
See [current validation and job status](V113_CLOUD_4000.md).

## Changes

- Store incoming electronic transitions as three continuous arrays: offsets,
  source indices and real-valued couplings. The Hamiltonian in this project
  has real star couplings. Each row is filled in ascending source order,
  preserving the previous scattered sum's accumulation order. The builder
  does not rely on Hermiticity to transpose the sparse graph.
- Build compact-mask neighbors directly using bit operations. The fermion
  sign counts occupied bath orbitals below the hopping orbital; set-based
  construction remains available for larger orbital counts.
- Accumulate H*psi per destination and immediately update the Runge-Kutta
  derivative. This removes the third full reference work vector. State and
  derivative must be distinct arrays; the existing integrator supplies them.
- Free graph-construction index maps before allocating propagation arrays.
- Correct path-local electronic-phase logging: no full phase table is allocated
  in this mode. Historical output printed its hypothetical full-basis size.

Physical parameters, Poisson random draws, reference selection, split steps,
measurement cadence and truncation levels are unchanged. The original tiny
amplitude cutoff is preserved, including at subnormal values.

## Local validation

```
python tests/test_reference_masks.py
python tests/test_reference_hopping.py
python tests/test_reference_rhs.py
```

652 determinant graph/order cases; 354 signed graphs / 44,148 directed edges
compared with explicit creation/annihilation on ordered sets and checked for
Hermiticity. 525 reference operator comparisons against tagged v1.12 cover
all five operator parts, zero/tiny/random complex amplitudes, accumulation
coefficients, nonsymmetric test graphs and uneven local oscillator slices
with halo rows. Maximum operator difference: 0.

The tests compile extracted production functions. The operator test also
exercises the actual new CSR builder; it requires the local v1.12 tag.
These local tests verify the changed graph and operator calculations; the
separate cloud report records full MPI regression and resource measurements.

For 715,136 states and 384 oscillator levels, the formula for global work
vectors falls from approximately 12.28 GiB (three arrays) to 8.18 GiB (two).
Halo buffers, reference metadata and sampling overhead are additional.
These are array-size calculations, NOT measured process memory or speedups.
