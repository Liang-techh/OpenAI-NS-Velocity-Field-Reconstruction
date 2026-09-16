# Amplitude log-scale error provenance

Status: **foundation and actual metadata integration implemented; focused coverage recorded**.
This artifact bounds only the scalar amplitude log-scale discrepancy. It does
not certify coefficient arithmetic, Bell derivatives, global `AxisSpace` norms,
or the paper-exact reconstruction.

## Pinned relation

The upstream pin is
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
`NavierStokes/NaturalAxisCoefficients.lean:47-49` defines

```text
realGradient(h,j,sigma,eta) = -L(h,eta)*H(h,j,eta)/(H(h,j,eta)^2+sigma^2),
```

`:413-414` defines `realPhase` as the oriented integral of that gradient, and
`:468-485` gives

```text
realAmplitude(eta) = exp(Lambda * realPhase(eta)) / C,
a'(eta) = Lambda * realGradient(eta) * a(eta).
```

The local signed-log amplitude state and selected symbolic `log_C` are recorded
in `references/AXIS_COEFFICIENT_AMPLITUDE_PROVENANCE.md`. The implemented
`src/openai_ns_reconstruction/axis_amplitude_log_scale.py` foundation adds the
rounding/error seam for its scalar log value. It has compiled and passed the
root arithmetic review. The exact foundation command
`python -m pytest -q tests/test_axis_amplitude_log_scale.py -W error` reports
`5 passed in 0.15 s`, including the multiplication-rounding case. The focused
integration command
`python -m pytest -q tests/test_axis_amplitude_log_scale_integration.py -W error`
reports `6 passed in 0.20 s`. The mixed-route command
`python -m pytest -q tests/test_axis_amplitude_log_scale_mixed.py -W error`
reports `9 passed in 0.22 s`, and the aggregate command
`python -m pytest -q tests/test_axis_amplitude_log_scale_aggregate.py -W error`
reports `4 passed in 0.18 s`. Together these are `24` distinct focused tests,
not full global acceptance.

## Exact interval-to-midpoint contract

Let `[L,U]` be the validated exact phase-to-log interval

```text
ell = Lambda * phase - log_C,
```

bound to the exact selected `h,j,sigma,eta`, finite selected `Lambda`, and
selected `log_C`. Let `ell0` be the actual Decimal midpoint used by the
amplitude state, and interpret that Decimal as the exact rational
`Fraction(ell0)`. The scalar log-scale discrepancy is enclosed by

```text
delta = [ L - Fraction(ell0), U - Fraction(ell0) ].
```

For an integer power `q >= 0`, let `s` be the actual Decimal scale used for
that power and interpret it as `Fraction(s)`. The corresponding discrepancy is

```text
delta_q = [ q*L - Fraction(s), q*U - Fraction(s) ].
```

This includes both phase-interval uncertainty and multiplication/Decimal
rounding in the actual `q`-power scale. The `q=0` scale is exactly zero. All
source inputs, phase bounds, midpoint conversion, and Decimal-to-Fraction
interpretations are exact in this contract.

Equal powers share the same base `delta` and therefore retain their common
correlation; they must not be assigned independent error intervals. Distinct
powers also depend on that same underlying `ell`, so their `delta_q` channels
are correlated through the common base interval. They remain separate because
each actual Decimal scale `s` contributes its own deterministic rounding offset;
the foundation does not yet perform a joint range evaluation or summation.

## Boundary

Integration of these intervals now reaches the actual amplitude, `x0`, all assigned
first-Picard nonlinear branches, the complete coefficientwise remainder/second-Picard path,
the formal triangular solver, and the wide profile metadata surfaces. The source is shared by
each signed-log channel and exact `Lambda`/eta identity checks fail closed before combination.
The focused metadata evidence totals `35` distinct tests: the prior `24`, plus `5` nonlinear
tests in `0.80 s`, `3` remainder/second-Picard tests in `0.38 s`, and `3` formal/wide-profile
tests in `0.22 s`. A separate existing-default-file batch returned `4 passed in 20.19 s` and
included one new nonzero-eta (`eta=-0.02`) metadata test, bringing the new metadata total to
`36` distinct tests. These tests do not certify exact Bell-factor/logarithm enclosures or full
coefficient arithmetic, and they do not establish SchedulePressure, parameter, global
weighted-space, or paper-exact claims.

The full suite immediately preceding this metadata work returned `1936 passed, 1 failed` in
`1723.45 s`; the sole fixture-argument issue was fixed and its affected test subsequently
passed in `0.20 s`. The full suite has not been rerun after that fix. The demo exited `0` with
all checks, and the direct strict-audit module exited `2` as expected according to its log; an
earlier console-wrapper exit `1` is historical. Relevant source compilation and review checks
passed, and the local checkpoint is being saved without an external push.
This contract does not bound Bell recurrence arithmetic, coefficient products,
pressure or SchedulePressure quadrature, selected-parameter admissibility,
global weighted-space membership, or any downstream profile or paper-exact
claim.
