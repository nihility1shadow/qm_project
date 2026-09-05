# v1.13 candidate: real incoming-edge CSR and two reference work vectors

## Status

Candidate only. Local graph and propagation-operator regressions pass.
The full Linux/MPI build, runtime/memory measurement and high-order long-time
accuracy tests have NOT run. Remote upload was blocked by automatic approval
review, which requires explicit authorization of the existing computation
server address. An approval question is pending with the user. Do not treat
this stage as proof of high fitting Q for 30 orbitals / 15 electrons.

The v1.12 compact-reference tag is the fully cloud-tested rollback point;
its completed 2500-a.u. accuracy limits are documented in V112_COMPACT_REFERENCE.md.
The launcher continues to default to v1.12. After building the candidate,
select it explicitly with BINARY=na_mpi_v113_csr.out.

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
This is not a substitute for compiling the entire simulation or testing real
MPI communication. It verifies the changed graph and operator calculations.

For 715,136 states and 384 oscillator levels, the formula for global work
vectors falls from approximately 12.28 GiB (three arrays) to 8.18 GiB (two).
Halo buffers, reference metadata and sampling overhead are additional.
These are array-size calculations, NOT measured process memory or speedups.

## Concrete validation plan after server approval

Use the existing project directory only; build a separate candidate executable.
First compare with the recorded same-seed v1.12 short outputs, then measure
D3 at 64 ranks and 128 ranks across two nodes. Example submissions:

```bash
sbatch --ntasks=4 --export=ALL,BINARY=na_mpi_v113_csr.out run_poisson_lowmem.slurm 200 1000 10 5 3 4096 1 111001
sbatch --ntasks=64 --export=ALL,BINARY=na_mpi_v113_csr.out run_poisson_lowmem.slurm 20 100 30 15 2 65536 1 111001
sbatch --ntasks=64 --mem=48G --export=ALL,BINARY=na_mpi_v113_csr.out run_poisson_lowmem.slurm 20 100 30 15 3 1000000 1 111001
sbatch --nodes=2 --ntasks=128 --ntasks-per-node=64 --mem=48G --export=ALL,BINARY=na_mpi_v113_csr.out run_poisson_lowmem.slurm 20 100 30 15 3 1000000 1 111001
```

Check correct observable files (ahm-sepmb-*.dat), particle conservation,
all finite values and same-seed differences. The original 64-rank D3 capacity
result was 369.30 s for 10 a.u. in v1.12; the candidate must be measured.
Only after these checks, select the faster rank count and run 2500 a.u.,
715,136 reference states, 384 oscillator levels with independent seeds,
within 48 hours per run. Assess each active orbital in every 100-a.u. window,
show all orbitals and distinguish seed repeatability from true QM fitting.
Retain a higher-Fock check and reference-order comparison. Do not smooth or
clip estimates to inflate fitting scores. No new cloud jobs have been submitted
at this candidate stage.
