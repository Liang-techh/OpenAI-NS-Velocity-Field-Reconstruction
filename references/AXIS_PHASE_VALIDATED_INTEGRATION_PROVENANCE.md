# Validated natural-axis phase integration provenance

Status: **provenance and conditional finite-cell theorem only**.  This artifact
does not claim a paper-exact amplitude, pressure, fixed point, or velocity.

## Pinned connection

The official pin is `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
`NaturalAxisData.lean:24-35` defines

```text
D = 1/2-h,  L(h,x)=1-2*h*x^2,
H(h,j,x)=D*x+(1-x^2)*(4*x+j).
```

`NaturalAxisCoefficients.lean:47-49` defines

```text
g(x) = realGradient(h,j,sigma,x) = -L(h,x)*H(h,j,x)/(H(h,j,x)^2+sigma^2).
```

Its `realPhase` definition at `:413-414` is
`eta * integral(t=0..1, g(t*eta))`, which equals the oriented integral
`integral(x=0..eta, g(x))` by the substitution `x=t*eta`, including negative
`eta`.  The amplitude definition and derivative theorem at `:468-485` give

```text
a(eta) = exp(Lambda*realPhase(eta))/C,
a'(eta) = Lambda*g(eta)*a(eta).
```

Here `C` is independent of `eta`; `realPhase(0)=0`, so the unresolved amplitude
offset is `log(a(0)) = -log(C)`.  The coefficient window is
`[-11/10,11/10]` (`NaturalAxisCoefficients.lean:23-27`), and `:65-78` supplies
the positive-real-axis `L` and denominator nonvanishing facts under the pinned
small-parameter and `sigma > 0` hypotheses.

The local symbols inspected are `natural_axis.py:real_gradient/real_phase`,
`axis_coefficient_amplitude.py:log_amplitude/jet_log`, and the rational input
and Bell provenance artifacts.  The phase integrand is rational for exact
rational `h,j,sigma`; `Lambda` and `C` are downstream scale inputs.

## Exact Fraction cell theorem

Let `[c-r,c+r]` be a cell with rational `c` and `r>0`, and put `rho=2*r`.
Expand the rational polynomials around `c`:

```text
P(z) = -L(h,z)*H(h,j,z) = sum P_k*(z-c)^k
Q(z) = H(h,j,z)^2+sigma^2 = sum Q_k*(z-c)^k.
```

All coefficients are exact Fractions and `Q_0 = H(c)^2+sigma^2 >= sigma^2`.
If

```text
theta = sum(k>=1) abs(Q_k)*rho^k / Q_0 < 1/2,
```

then `|Q(z)| >= Q_0*(1-theta) > Q_0/2` on `|z-c| <= rho`.  With the exact
rational bound
`M = 2*sum(k>=0) abs(P_k)*rho^k/Q_0`, Cauchy's estimate gives Taylor
coefficients `b_n` of `g=P/Q` satisfying `abs(b_n) <= M/rho^n`.  They are
computed without division by a sampled value by

```text
b_0 = P_0/Q_0,
b_n = (P_n - sum(k=1..n) Q_k*b_(n-k))/Q_0.
```

The exact degree-`N` cell integral is

```text
I_N = sum(n=0..N) b_n * mu_n,
mu_n = 2*r^(n+1)/(n+1) for even n, and 0 for odd n.
```

The Cauchy tail on the cell is bounded by

```text
E_N = 2*r*M*2^(-N-1)/(1/2).
```

Therefore `[I_N-E_N,I_N+E_N]` encloses the cell integral.  Every endpoint,
coefficient, partial sum, and error is rational; convert the final interval to
Decimal with outward `ROUND_FLOOR`/`ROUND_CEILING` only after the exact sum.

Partition `[0,eta]` into rational symmetric cells for `eta>0`; for `eta<0`,
partition `[eta,0]` and negate the summed interval.  `eta=0` is exactly zero.
Allocate positive rational budgets `tau_i` with `sum tau_i <= tau`, and raise
each cell's order until `E_N <= tau_i`.

Refining cells alone at fixed `N` is not a termination proof: since
`E_N/(2*r)=M*2^(-N)`, the sum over a fixed-length partition need not decrease.
Order growth is required, and a finite resource cap must fail closed.

For fixed positive rational `sigma` on a compact real interval `K`, finite
termination follows conditionally.  Since `Q(x)>=sigma^2`, and `Q` is a fixed
polynomial, there is a finite uniform positive rational majorant `B` (take
`B=1` if the nonconstant coefficient sum vanishes) with
`sum(k>=1) abs(Q_k(c))*rho^k <= B*rho` for every `c in K` and `rho<=1`.
Choosing cells eventually with `rho < sigma^2/(2*B)` makes `theta<1/2`; exact
bisection therefore reaches the gate after finitely many splits.  On each
accepted cell, `M` and `r` are finite and `E_N` tends to zero as `N` increases,
so the allocated tolerance is reached after finitely many order increments.
Small `sigma` can force very small cells near the cubic `H` root, but does not
create a real pole.  The proof assumes exact rational inputs and exact gate
comparisons; interval parameter inputs require an interval version of the same
inequalities.

## Dyadic outward accumulation contract

The exact cell theorem above can be accumulated without constructing a large
LCM across all cell denominators.  For a positive-oriented cell, let

```text
b = tau * cell_length / total_length
```

be its exact rational budget, and accept the cell only after the exact local
Cauchy error satisfies `E <= b/2`.  Choose an exact dyadic quantum
`d = 2^q <= b/4`; `q` may be positive for a large budget or negative for a
small one.  Round the exact cell interval outward to

```text
lo = floor((I-E)/d) * d,
hi = ceil((I+E)/d) * d.
```

The true cell integral remains in `[lo,hi]`.  Since the exact interval radius
is `E`, its rounded radius obeys

```text
(hi-lo)/2 <= E+d <= b/2+b/4 = 3*b/4.
```

Summing cell intervals therefore gives total radius at most
`3/4 * sum(b) <= 3*tau/4`, while every accumulated endpoint is dyadic.  Sums
of dyadics have a common power-of-two denominator, so the accumulator avoids
the large general-rational LCM.  This is a rounding enclosure, not a claim
that the rounded midpoint is the exact Taylor sum; all local `I` and `E`
calculations remain exact/certified before projection.

The implementation keeps the public exact-Fraction accumulation mode for
compatibility and selects this dyadic outward mode for the default amplitude
path.  This does not weaken the requested tolerance or change any
`h,j,sigma,eta` input.  A
cell that fails the denominator gate, `E <= b/2`, or any finite order/cell/depth
cap fails closed before rounding.  For negative `eta`, accumulate sorted
positive cells and negate the final interval as `[-upper,-lower]`; `eta=0`
returns `[0,0]`.  Exact integer floor/ceiling handles negative `I` and dyadic
quanta outside binary64 range.

## Bounded performance evidence

A bounded diagnostic at the selected `eta=-1/50` completed the exact rational
phase path in `10.850 s`, with `225` accepted cells and maximum order `64`.
The affine log-interval transport took `0.558 s`.  The Decimal midpoint path
measured `11.307 s`, with stack evidence in Decimal/Fraction conversion.  The
current default validated binary64-input diagnostic began at `23.019 s` and hit
the `45 s` bound before completing; no precise runtime is claimed.

The completed slow2 validation recorded
`python -m pytest -q tests/test_axis_phase_integral.py -W error` as `7 passed
in 0.20 s`.  The combined command
`python -m pytest -q tests/test_axis_amplitude_validated_default.py
tests/test_axis_coefficient_amplitude.py tests/test_axis_phase_log_enclosure.py
-W error` returned `13 passed in 43.71 s` under the hard `120 s` cap.  This
resolves the default-regression performance blocker for that bounded command.
The earlier approximately `9` minute interruption remains historical and had
no pytest summary; no full-suite or paper-exact claim follows.

## Boundary

This validates only the chosen rational kernel integral.  It does not enclose
the theorem-selected real `h,j,sigma` if they are merely binary64 approximations,
the wide `Lambda`, `C` or `log(C)` offset, the SchedulePressure datum, amplitude
magnitude, accumulated coefficient arithmetic, global `AxisSpace` norms, or
the final reconstruction.  The exact-Fraction public integral mode remains a
compatibility path; the default amplitude mode is intended to use the dyadic
outward accumulator once its implementation baseline lands.  The selected
finite `Lambda` and symbolic `C` affine bridge remains conditional on those
scalar representations and does not itself certify their theorem-side
parameter selection.  The next Stage 1 boundary is propagating the phase
interval uncertainty through the genuine wide natural-remainder/profile chain
and independently bounding coefficient arithmetic errors; global weighted
`AxisSpace` certification remains downstream.
