# 8M trajectory parameter search protocol

- Paths per independent run: 8,000,000
- Independent repeats per parameter: 3
- Orbitals/electrons: 10/5
- Time range: 0-150 a.u. with dt=0.5 a.u.
- Real measurement stride: 1 output step (0.5 a.u.)
- Ranking interval: 0-130 a.u.
- Segment width: 10 a.u.
- Pass condition: minimum active-orbital Q > 4 in every segment ending at or before 130 a.u.
- Selection: retain ceil(30% of candidates) after every round
- Tie-breakers: wider physically meaningful wc, lower inactive-orbital noise, then higher Q in 120-130 a.u.
- All evaluated and rejected rows remain in the cumulative catalog.
