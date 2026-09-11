# Outgoing tail / clock-weight provenance

Status: **formal-structure subconstruction**, not `paper-exact`.

Pinned source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Mapped modules:

- `NavierStokes/OutgoingSchedule.lean`: `Parameters`, `logAmplitude`, `radialAmplitude`, `shape`, `angular`.
- `NavierStokes/OutgoingTail.lean`: `TailData`, `tailShape`, `tailDebt`, `releaseSlope`, `releaseLag`, `decayHold`, `flattenFactor`, `flattened`, `releaseAdjustment`, `tailStart`, `finalAngular`, `powerConstant`.
- `NavierStokes/SchedulePressure.lean`: `clockWeight`, `shapeExponent` and the exact factorization used later for `axisPressure`.

## Constructive choice

Pinned Lean defines `OutgoingTail.stepBound` by `Classical.choose`. The repository previously proved that `S = 32` is an admissible global upper bound for the pinned smooth-step derivative. This module consistently instantiates

`flattenLength = 10 * (32 + 1) * log(2) + 1`

and

`tailCoefficient = exp(-5) / (16 * (32 + 1))`.

This is a theorem-admissible constructive realization of the same printed formulas. It is **not** a claim that Lean's opaque chosen value is definitionally equal to 32.

## What is executable now

`src/openai_ns_reconstruction/outgoing_tail.py` materializes the scalar outgoing schedule through the complete `finalAngular(y, eta)` and `clockWeight(y) = finalAngular(y,0)^2` for any finite parameters satisfying the exact `OutgoingSchedule.Parameters` / `TailData` inequalities. Long plateaus are reduced algebraically to the pinned smooth-step primitive, so numerical quadrature is confined to compact transition intervals rather than being used as a fit over the whole schedule.

The release integrating-factor computation is independently regression-checked against `scipy.integrate.solve_ivp`; the final angular evaluator is also checked against the exact ideal-prefix identity, eta-independence after flattening, and the eventual power-tail identity.

## Remaining boundary

This does **not** yet materialize `SchedulePressure.axisPressure = -(1/2) integral_R finalAngular(y,eta)^2 dy`; therefore it does not yet certify the low-|Z| uniform positive `H^2` margin, choose the resulting natural-axis cutoff sigma, construct Lambda/C, or instantiate the coefficient-space fixed-point fields. No leading profile is promoted to `paper-exact` by this increment.
