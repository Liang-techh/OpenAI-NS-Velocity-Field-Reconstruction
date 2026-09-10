# Numerical implementation and limits

## 1. Concentration scale and time representation

The original bisection stopped when the bracket width was at most
`rtol * max(1, q)`. For small `q`, this was a unit-scale absolute tolerance,
not a relative tolerance. The stored regression case uses `h=0.005`,
`eta=0.6`, and the representable `tau=1-(1-1e-14)`. Its old relative error
was about 1.46044; the new error is approximately 3.05e-14.
The exact input/output numbers are in `reports/benchmark.json`.

For `k=1-2h`, choose `scale=max(tau, |z|^(2/k))` without first forming the
potentially overflowing power. With `q=scale*exp(y)`, solve for `y>=0`:

```text
1 - exp(log_gamma-k*y) - exp(log_a-y) = 0,
log_gamma = k*(log_axis-log_scale),
log_a = log(tau)-log_scale.
```

The implementation uses `expm1`, brackets on the physical branch, and stops
on a relative logarithmic width. It reports failed bracketing, overflow,
invalid arguments and nonconvergence. Float64 rounding of the inputs and
transcendental functions is not enclosed by interval arithmetic.

`*_from_tau` avoids the separate problem of `1-tau` rounding to 1. Preserve
`d=tau/q` instead of subtracting `eta**2` from 1 near the endpoints. A
rounded `eta` very close to 1 can coexist with a much smaller positive `d`.
`d` is the stable value used in the differential identities.

## 2. Radial averages and axis behavior

Using a unit-interval substitution,

```text
A_X(U) = integral_0^1 U(X*s,eta) ds.
```

This avoids small-integral/small-X division. A cached 32-node Gauss–Legendre
rule replaces the old 801-point trapezoidal rule. Optional `average_U` and
`average_dU_deta` callbacks allow an independently known antiderivative.
**32 nodes are not sufficient for arbitrary high-frequency profiles.**
Callers must compare increased quadrature orders and check the independent
identities needed by the construction.

For a smooth axisymmetric Cartesian field, `E=sqrt(2X)*F` with smooth `F`.
An optional `F` callback supplies its axis limit and gives `Pi_X=F^2`.
Callers supplying both `E` and `F` must make them consistent; the interface
cannot prove equality of arbitrary callbacks. The Gaussian fixture does
so analytically. A missing pressure profile is an error, not an implicit
choice of zero pressure. Nonzero radii below 1e-15 are no longer erased.

## 3. Analytic background curl, including cutoff terms

For one coefficient with `lambda=2*n*h`:

```text
A_n = 0.5*q^(-A+lambda)*A_X(U_n)*(-y,x,0),
v_lambda = [2*eta*U_n - 2*(D+lambda)*eta*A_X(U_n)
            - d*d_eta(A_X(U_n))] / L.
```

Writing `w=chi(a*q)` and `wq=a*chi'(a*q)`, the radial velocity divided by
radius is

```text
u_r/r = 0.5*w*q^(lambda-1)*v_lambda
        - wq*eta*q^lambda*A_X(U_n)/L.
```

The second term is essential in the transition region. It comes from
`grad(w) cross A_n`; multiplying an already computed poloidal velocity by
`w` would omit it. Tests compare this analytic evaluator to an independent
finite-difference curl of the potential, including n=1,2,5 and cutoff
transition points. Zero-order assembly is checked against leading velocity.

These are finite sums of caller-supplied profiles, **not** a solver for the
Section 5 recursion or a validated infinite series. The C-infinity cutoff
is an explicit experimental bump, not a claim to have selected the paper's
recursive admissible scales. The base n=0 term is not cut off here.

## 4. Appendix A moment primitives

`power_moment_matrix` forms a finite moment matrix from distinct powers and
normalized smooth bumps on ordered disjoint positive intervals. It uses
numerical quadrature; a determinant or singular value is not a rigorous
uniform invertibility proof.

`solve_quadratic_moments` implements the zero-initialized iteration described
in Lemma A.2:

```text
B*c + Q(c,c) = d,
c_next = solve(B, d-Q(c,c)).
```

It checks the sufficient condition `8*beta^2*kappa*||d|| <= 1`, where beta
is estimated from singular values and the Frobenius norm bounds the
bilinear norm kappa. For `Q=0`, the smallness condition is unnecessary.
The result includes the residual, iterations, smallness estimate and ball
radius. Invalid shapes, ill-conditioning, a failed smallness screen, and
nonconvergence are not silently accepted.

This implementation handles **one parameter value at a time**. It does not
provide the actual matching data, eta-uniform bounds, parameter derivatives,
profile cutoffs, shear modification, or cone conditions. All bounds are
floating-point estimates rather than formally certified enclosures.

## 5. Localization and independent diagnostics

A localized potential gives

```text
curl(c*A) = c*curl(A) + grad(c) cross A.
```

An analytic curl and cutoff gradient can avoid nested finite differences.
The poloidal curl is divergence-free as a differential identity. The direct
term `c*B*e_theta` additionally requires theta-independence of `c*B`.
A negative regression test shows that a generic nonaxisymmetric cutoff
breaks this property. A generic callback is not declared divergence-free.

The momentum convention is

```text
u_t + (u.grad)u - viscosity*Delta(u) + grad(p) = f.
```

Tests include a polynomial manufactured solution for which all four terms
are nonzero and the force is independently given, a wrong-force rejection,
and an independently specified unforced Taylor–Green solution. Refining
space/time steps by two gives approximately fourfold smaller residuals.
Computing `f` with the same residual routine and subtracting it would only
check a tautology; it does not certify smooth forcing through t=1.

Time stencils accept an optional open domain such as `(None,1.0)` and shrink
without crossing its boundary. Supply analytic derivatives when a valid
Float64 stencil cannot be represented. Very small steps can amplify
roundoff, especially when numerically differentiating a numerical curl.

`kinetic_energy_axisymmetric` includes the `2*pi*r` cylindrical measure, but
integrates only a finite cylinder at one time. It is not a whole-space or
uniform-in-time energy estimate. The exported toy slice is neither a CFD
time evolution nor a reconstruction of the paper's final compact field.
