# Conditional amplitude derivative enclosures

Status: implemented and reviewed for the chosen rational kernel and selected
finite scalars; focused acceptance passed. No global reconstruction claim.

This extends the pinned amplitude relation documented in
`AXIS_COEFFICIENT_AMPLITUDE_PROVENANCE.md`, for
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
It does not establish admissibility of the selected SchedulePressure scalars,
global coefficient-space bounds, or a paper-exact velocity field.

## Derivative identity

Write `ell = Lambda * realPhase - log_C` and `g = realPhase'`.
For `q = 1, 2`, the full parameter derivative is

    d_eta^m exp(q ell) = exp(q ell) B[q,m].

The exact rational Bell recurrence already implemented by
`RationalAxisCoefficientData.amplitude_power_bell_fraction` uses the exact
rational kernel derivatives and the exact rational interpretation of the
selected finite Decimal Lambda. In particular,

    B[q,0] = 1
    B[q,1] = q Lambda g
    B[q,2] = q Lambda g' + (q Lambda g)^2.

There is no assumption that `g(0)` vanishes. Nonzero j can make it nonzero.
The sign is the exact sign of B; B = 0 yields an exact zero derivative.

For a nonzero B, an enclosure `[L,U]` for ell and `[l,u]` for
`log(abs(B))` imply the full derivative log-magnitude enclosure

    [q L + l, q U + u].

All four endpoints are rational, so this addition loses no small correction
beside a very large amplitude log scale. The factor is not first rounded to
Decimal before its sign or logarithm is enclosed.

## Rational logarithm proof

For x > 0 choose an integer k exactly such that x = 2^k m and 1 <= m < 2.
Then log(x) = k log(2) + log(m). For t = (m-1)/(m+1),

    log(m) = 2 sum(i=0..N-1, t^(2i+1)/(2i+1)) + R_N,
    0 <= R_N <= 2 t^(2N+1) / ((2N+1)(1-t^2)).

The same formula for log(2) uses t = 1/3. A negative k reverses the
endpoints of the log(2) interval when it is multiplied by k.

To avoid repeatedly raising the original Bell fraction's denominator to high
powers, enclose t on a dyadic grid of step d <= min(1/16, b/16).
The rounded bounds obey 0 <= t_lower <= t <= t_upper <= 1/2.
On this interval the derivative of `2 atanh(t)` is at most 8/3. Thus the
uncertainty from rounding t is at most b/6. Choose N so the upper-endpoint
tail is at most b/2. The lower partial sum and upper partial sum plus tail
then enclose the logarithm. Rounding these two endpoints outward to the same
grid adds at most 2d <= b/8 to their width. Total width is at most

    b/6 + b/2 + b/8 = 19b/24 < b.

For requested total log width epsilon, use b = epsilon/(abs(k)+1).
Combining the two component intervals gives width at most epsilon. Exact
comparisons check the budget, and a finite term cap raises an error if the
requested bound is not achieved. No numerical convergence guess or ordinary
rounded logarithm serves as a certificate.

## Limits of this result

The logarithm tolerance controls the Bell-factor logarithm only. Total
derivative uncertainty also includes q times the phase-log interval width;
large Lambda can make this contribution enormous. This is a bound for the
chosen rational kernel and chosen finite scalars, not a new proof of their
selection. General pressure/coefficient products and sums still need their own
arithmetic enclosures. Existing point-valued jet APIs retain their separate
numerical status.

## Measured acceptance

`ActualScheduleAmplitudeLogState.derivative_log_enclosure(power, order, eta)`
uses the exact Bell factor and the bounded rational logarithm above. The
standalone result constructor remains conditional on its supplied factor and
log interval; the actual amplitude factory computes those inputs.

`python -m pytest -q tests/test_axis_amplitude_derivative_enclosure.py -W error`
returned **5 passed in 2.32 s**. Coverage includes independent high-precision
logarithm references, reciprocal orientation, width and finite-cap checks,
the explicit first two derivative identities for actual data at eta zero,
and exact-zero representation. Source compilation passed. These focused
results do not replace a full current-tree regression or global proof.
