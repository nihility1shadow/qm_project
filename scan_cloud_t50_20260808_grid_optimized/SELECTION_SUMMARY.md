# 50t optimized parameter selection

## Fixed scope

- Norb=10, Nel=5, dt=0.5 a.u., t=0..50 a.u. (101 points).
- 2,000,000 trajectories per independent run, 16 MPI ranks.
- EF=-4.5 eV; E0=EF-wc/2; E1=EF+wc/2.
- O 2p delE=-17.194874968839816 eV.
- Search limits respected: 0 < wc <= 10 eV; 3e-5 <= eta <= 2.2e-4.

## Selection

Primary: wc=6.5 eV, eta=7.5e-05.
Three-repeat Q segments: 28.225, 7.966, 3.937, 2.290, 2.259.
Qconv=2.259; repeat correlation min/mean=0.876/0.889.
Conservative N4000 estimate: Q=3 -> 2.82e+08; Q=5 -> 7.84e+08 trajectories.

## Interpretation

Q is RMS repeat-mean orbital signal divided by RMS between-run Monte Carlo noise in each 10 a.u. segment.
The 4000t estimate assumes Monte Carlo noise scales as N^(-1/2) and accumulated difficulty as sqrt(t). It is a planning estimate, not a 4000t validation.
Primary ranking uses the minimum segment Q, so a candidate cannot win only because its final segment happened to be quiet.

## Files

- optimized-50t-comparison.png: segmented Q, trajectory estimate, best orbital changes, and signal/noise.
- best-all-orbitals-three-repeats.png: all ten orbitals for the three independent best runs.
- segmented_convergence.csv/json: complete numeric results for the four retained candidates.
- Candidate subdirectories: retained raw SepMB output and run metadata only.
