# Actual SchedulePressure theorem-selected realAmplitude log state

Status: **formal-structure only**. This artifact does not make the Stage-1
leading profile paper-exact and does not change
`paper_exact_velocity_available=false`.

## Pinned source

This implementation follows
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
especially `NavierStokes/NaturalAxisCoefficients.lean`:

- `realGradient(h,j,sigma,eta) = -L H / (H^2 + sigma^2)`;
- `realPhase` is the integral primitive of `realGradient` with
  `realPhase(0)=0`;
- `realAmplitude = exp(Lambda * realPhase) / C`.

The actual SchedulePressure chain in this repository already selects the
theorem-side `Lambda` and symbolic normalization threshold
`C = exp(C_exponent)` in wide `Decimal` arithmetic. This increment binds those
selections to the actual schedule and represents the resulting amplitude
without narrowing either scalar to binary64.

## Executable representation

`axis_coefficient_amplitude.py` introduces
`ActualScheduleAmplitudeLogState`. Its factory accepts only the actual
`TailData` and `j`; it internally reconstructs the landed
`ActualScheduleReferenceAxisState`, fixed `AxisData`, and wide scale chain.
There is no caller override for pressure data, `sigma`, `epsilon`, `Lambda`,
`C`, or coefficient values.

For radial degree zero it stores each parameter jet as a sign and logarithmic
absolute magnitude. The zeroth jet uses

`log a(eta) = Lambda * realPhase(eta) - C_exponent`.

Higher eta derivatives use the exact structural identity

`a' = Lambda * realGradient * a`

and the complete Bell recurrence. Radial rows `n>0` are exactly zero, matching
the fact that `realAmplitude` depends only on eta.

A `binary64_state()` adapter is provided only as a fail-closed bridge to the
existing coefficient backend. Any mathematically nonzero jet below binary64's
minimum subnormal magnitude, or above its maximum finite magnitude, raises
`ArithmeticError`; it is never silently replaced by zero or infinity.

## Validation and boundary

Regression coverage checks that the state is bound to the actual wide
`Lambda/C` selection, verifies `log a(0) = -log C`, checks the first derivative
identity in log magnitude, verifies exact radial-degree-zero behavior, and
requires binary64 range loss to fail closed on the current conservative scale.

The remaining boundary is substantive:

1. `realPhase` values are still evaluated by the landed floating-point
   quadrature, so this is not an interval or Lean certificate.
2. The existing `AxisCoefficientJetState`/`naturalRemainder` backend is
   binary64-valued. The current theorem-selected amplitude and
   `t=1/Lambda` can lie outside that dynamic range. A wide/log coefficient
   arithmetic bridge is therefore required before evaluating the genuine
   `naturalRemainder(x0)` without losing nonzero terms.
3. The Picard fixed point, derived average/pressure fields,
   `NaturalProfileAssembly`, and final support/moment/matching/cone checks
   remain unresolved.
