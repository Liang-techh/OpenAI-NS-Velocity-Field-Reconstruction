# Rational amplitude power-Bell provenance

Status: **exact finite normalized-derivative algebra only**.  This note does
not certify the amplitude magnitude, pressure, a global norm, or velocity.

## Pinned identity

At commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
`NaturalAxisCoefficients.lean:413-414` defines `realPhase` as the integral of
`realGradient`, while `:468-485` defines

```text
a(eta) = exp(Lambda * realPhase(eta)) / C,
a'(eta) = Lambda * realGradient(eta) * a(eta).
```

`C` is a constant in `eta`; it is the normalization denominator, not a
parameter-dependent factor.  The cached implementation is
`axis_coefficient_amplitude.py:ActualScheduleAmplitudeLogState._bell_factor`.

## Exact normalized Bell factors

For an integer `p >= 0`, put `f = log(a)` and
`q_r = p * Lambda * realGradient^(r-1)(eta)` for `r >= 1`.  Then

```text
(a^p)^(m) / a^p = B_m(q_1, ..., q_m),
B_0 = 1,
B_(m+1) = sum(k=0..m) binom(m,k) q_(k+1) B_(m-k).
```

The rational input provider from `AXIS_COEFFICIENT_RATIONAL_DATA_PROVENANCE.md`
gives every `realGradient^(r)(eta)` exactly: form the Taylor series of
`-L(eta+t) H(eta+t)/(H(eta+t)^2+sigma^2)` with Fraction products and
reciprocals, then multiply coefficient `r` by `r!`.  For exact rational choices
of `h,j,sigma,Lambda`, every `q_r` and `B_m` is rational.  The factor `C`
cancels from this normalized ratio, as does the phase value; the phase integral
and `C` still determine the unnormalized magnitude `a^p`.

For each exact rational `B_m`, a directed Decimal enclosure is obtained by
outward `ROUND_FLOOR` and `ROUND_CEILING` conversion at the selected precision.
This closes input and Bell-factor roundoff only.  It does not enclose
`realPhase`, `log(a)`, wide amplitude scale, or subsequent coefficient products.

## Existing comparison and boundary

Before this change, the legacy recurrence used binary64 gradient jets.  The
current `axis_coefficient_amplitude.py`, wide-source, and formal-solver paths
use the rational Bell provider; the phase value remains supplied by numerical
quadrature.  The old result was a point calculation, not an exact or interval
Bell value.

The next phase primitive is a validated enclosure for
`realPhase(eta) = integral_0^eta realGradient(x) dx`.  The integrand is rational
for rational `h,j,sigma`, but a negative path can cross the narrow cubic `H`
root near `(-j/4,-j/5)`; fixed 4001-point trapezoids are not an error bound.
Adaptive outward-rounded rational quadrature, split at or around that root,
would certify the phase without changing the pinned formula.  No global
parameter, pressure, or paper-exact claim follows from this Bell layer.
