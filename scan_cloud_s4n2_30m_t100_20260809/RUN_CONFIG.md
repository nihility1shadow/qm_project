# Two-run 30-million-trajectory validation

- `Norb=4`
- `Nel=2`
- `ntraj=30000000` per independent run
- independent runs: 2
- `wc=0.25 eV`
- `eta=6e-5`
- `AHM_DELE_EV=-2.5`
- `AHM_NSTEP=200`
- `SEP_MB_TMAX=100`
- `SEP_MB_MEASURE_STRIDE=1`
- output interval: 0.5 a.u.
- MPI ranks per run: 64

The one-million-trajectory baseline had strict `Q(0..100)=2.1841` over six
independent repeats. Square-root scaling predicts approximately `Q=11.96` for
30 million trajectories.
