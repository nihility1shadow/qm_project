# Poisson sampling noise and optimization direction

## 1. Scalar Poissonization

For a constant Hamiltonian `H = hbar * omega`, introduce a Poisson event count
`N_t ~ Poisson(lambda * t)` and the unbiased stochastic propagator

`U_st(t) = exp(lambda*t) * (-i*omega/lambda)^(N_t)`.

Using `E[z^N] = exp(lambda*t*(z-1))` gives

`E[U_st(t)] = exp(-i*omega*t)`

and

`E[|U_st(t)|^2] = exp[t*(lambda + omega^2/lambda)]`.

Therefore

`Var[U_st] = exp[t*(lambda + omega^2/lambda)] - 1`.

The exponent is minimized by

`lambda_opt = |omega|`,

but even at the optimum the variance is

`Var_min = exp(2*|omega|*t) - 1`.

Rate optimization minimizes the exponential growth coefficient; it does not
remove exponential long-time variance.

## 2. Matrix and many-body form

For `H = H0 + V`, Poissonize the interaction-picture perturbation `V_I(t)`.
A trajectory with event times `tau_1,...,tau_K` carries a complex product of
interaction matrix elements divided by its sampling rates. A norm bound for a
constant total rate `Lambda` has the same structure as the scalar result:

`second moment <= exp[t*(Lambda + G^2/Lambda)]`,

where `G` is a characteristic norm or spectral radius of `V/hbar`. This bound is
minimized at `Lambda = G`, with irreducible exponent `2*G*t`.

For the homogeneous central-orbital many-body graph in this project,

`c = sqrt(eta/(Norb-1))`,

`Nvac = Norb - Nel`,

and the two electronic classes have `Nel` and `Nvac` possible jumps. Their
balanced spectral scale is

`G = sqrt(Nel*Nvac) * |c|`.

The code uses exactly

`rate = G * dt`.

The jump compensation factors `sqrt(Nvac/Nel)` and `sqrt(Nel/Nvac)` keep a
forward-and-return pair balanced. Thus changing only the homogeneous Poisson rate
cannot provide an order-of-magnitude convergence improvement.

## 3. Why the observable is noisier than an amplitude

Each occupation is a forward-backward projection. A sampled contribution has the
form

`X = Re[w_for * conj(w_back) * overlap]`.

The phase contains `(-i)^(K_for-K_back)`, electronic dynamical phases, coherent
state phases, and a fermionic permutation sign. Individual `|X|` values can be
large while positive and negative contributions nearly cancel in the mean. This
is the dynamical sign/phase problem.

The reported observable is also a ratio,

`A_hat = X_bar / Z_bar`,

where `Z_bar` is the sampled trace normalization. To first order,

`Var[A_hat] approximately Var[X_bar]/Z^2 + A^2*Var[Z_bar]/Z^2`

`                       - 2*A*Cov[X_bar,Z_bar]/Z^2`.

Cancellation noise in both numerator and denominator is therefore inherited by
the normalized occupation.

For a determinant at graph distance `d`, at least `2*d` jumps, or `2*d+1` when
the electronic class changes, are required. The relevant Poisson probability is

`P(K=k) = exp(-G*t) * (G*t)^k/k!`.

Rare high-order sectors have small physical means but large importance weights.
This is why weakly occupied orbitals can have the worst relative noise.

## 4. Numerical scales in the two configurations

For the four-orbital scan (`Norb=4`, `Nel=2`, `eta=6e-5`):

`c = sqrt(6e-5/3) = 0.0044721 Ha`,

`G = sqrt(2*2)*c = 0.0089443 Ha`.

At `dt=0.5`, the event probability is `G*dt=0.0044721`, and the expected event
count at `t=100` is `G*t=0.8944`.

For the restored oxygen default (`Norb=10`, `Nel=5`, `eta=0.01`):

`c = sqrt(0.01/9) = 0.0333333 Ha`,

`G = sqrt(5*5)*c = 0.166667 Ha`.

At `dt=0.5`, the event probability is `0.0833333`, and the expected event count at
`t=100` is `16.6667`. The number of phase-cancelling histories is consequently
much larger than in the four-orbital weak-coupling scan.

The oxygen atomic `2p` energy is

`delE = -0.6319 Ha = -17.1948749688 eV`.

With the default `wc=10 eV`, bath energies lie approximately between `-9.5 eV`
and `+0.5 eV`. The central orbital is strongly off resonance. First-order transfer
has the scale

`|A_j(t)| approximately 2*|c_j|*|sin(Delta_j*t/2)|/|Delta_j|`.

Large detuning suppresses the physical occupation change, while the stochastic
weight variance remains. Absolute noise may decrease, but relative convergence
of a tiny signal can still be difficult.

## 5. What parameter tuning can and cannot do

- `wc` mainly changes bath energies and phase cancellation. It does not change
  homogeneous `c` in `diseven`, so it cannot directly remove weight growth.
- `eta` changes `c` and therefore `G`, but it also changes the physical
  Hamiltonian and the true transition amplitude. Lowering `eta` is not a pure
  numerical optimization.
- `Ntraj` only gives the ordinary Monte Carlo improvement
  `standard error proportional to 1/sqrt(Ntraj)`.

## 6. Recommended optimization order

1. Keep physical `wc`, `eta`, and `delE` fixed.
2. Split `H = Href + DeltaH` so that as much dynamics as possible is propagated
   analytically or deterministically, then Poissonize only `DeltaH`. The relevant
   variance scale becomes `||DeltaH||`, not `||H||`.
3. Use guided transition rates
   `r(x->y) proportional to |V_yx|*h(y)/h(x)`, where `h` approximates the magnitude
   of the remaining contribution. This is a Doob-transform importance sampler.
4. Stratify by jump order, endpoint sector, and Johnson/Kondo distance so rare
   orbitals receive controlled sample counts.
5. Couple forward and backward paths instead of sampling them almost
   independently; correlated phases reduce variance of their product.
6. For long times, use fixed-length windows with resampling, inchworm stitching,
   or sequential Monte Carlo. Direct Poissonization retains exponential
   long-time variance even at its optimal rate.

For the requested 4000 a.u. regime, the windowed/inchworm direction is essential.
Increasing the number of direct trajectories alone cannot reliably defeat an
exponential-in-time second moment.
