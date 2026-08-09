# Norb=4, Nel=2, 1M trajectory convergence search

- `Norb = 4`, `Nel = 2`
- `ntraj = 1,000,000` per independent run
- `t = 0..100 a.u.` with real measurements every `0.5 a.u.`
- `AHM_DELE_EV = -2.5 eV`
- three independent repeats per exploratory candidate
- Q windows are 10 a.u. wide
- strict score is the minimum active-orbital Q over all windows through 100 a.u.
- target is strict `Q > 10`
- active orbitals are classified from `t = 0..40 a.u.` at 5% of peak RMS

