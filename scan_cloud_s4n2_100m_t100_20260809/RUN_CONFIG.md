# One-run 100-million-trajectory validation

- `Norb=4`
- `Nel=2`
- `ntraj=100000000`
- independent runs: 1
- `wc=0.25 eV`
- `eta=6e-5`
- `AHM_DELE_EV=-2.5`
- `AHM_NSTEP=200`
- `SEP_MB_TMAX=100`
- `SEP_MB_MEASURE_STRIDE=1`
- output interval: 0.5 a.u.
- MPI ranks: 64

This run is compared against the mean of two independent 30-million-trajectory
runs. A single run cannot independently estimate repeat-to-repeat strict Q.
