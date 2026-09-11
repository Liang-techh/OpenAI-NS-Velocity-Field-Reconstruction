# SchedulePressure axis-pressure provenance

Truth status: **formal-structure / numerical evaluator of the actual constructed datum**. This file records a real implementation of the pressure integral attached to the landed outgoing schedule. It is not a proof of the later uniform low-|Z| margin, and it does not make the leading profile paper-exact.

## Pinned sources

Official Lean repository: `openai/NavierStokesAndEuler` at commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Primary module: `NavierStokes/SchedulePressure.lean`.

The implementation uses the following pinned statements/definitions:

- `SchedulePressure.axisPressure d eta := -(1/2) * ∫ y, OutgoingTail.finalAngular d (y,eta)^2`.
- `SchedulePressure.angular_square_factorization`, which identifies the same integrand with `clockWeight d y * PressureDatum.kernel (shapeExponent d y) eta`.
- `SchedulePressure.clockWeight_ideal` and `shapeExponent_ideal` on `y <= 0`.
- `SchedulePressure.axisPressure_hasDerivAt`, whose derivative contains `∫ shapeExponent * finalAngular^2`.
- `OutgoingTail.finalAngular_eventual_power`, which gives the exact exponential tail for `y >= tailEnd`.
- `OutgoingSchedule.radialAmplitude_hold` and the printed release/tail plateau formulas, used only to integrate intervals whose logarithmic slope is exactly constant.

`NavierStokes/PressureDatum.lean` supplies the real kernel

`kernel(a,eta) = exp(-2*a*log(1+eta^2))`.

## Executable artifact

`src/openai_ns_reconstruction/schedule_axis_pressure.py` evaluates the actual all-real-line datum for an executable `TailData`. It does not substitute the earlier `IdealPrefixPressureDatum` positive-clock continuation.

The improper integral is decomposed into mathematically exact pieces whenever the pinned construction is on a constant-slope plateau:

- `(-∞,0]`: exact ideal-prefix mass `5 P^2 (1+eta^2)^-2`;
- constant radial-amplitude-slope plateaus: exact exponential antiderivatives;
- `[tailEnd,∞)`: exact consequence of `finalAngular_eventual_power`;
- the six genuine smooth transition intervals: cached Gauss-Legendre quadrature of the factorized paper integrand.

The module also evaluates the derivative formula from `axisPressure_hasDerivAt`. The exponent-weighted integral simplifies further because `shapeExponent=1` before `endpoint` and `shapeExponent=0` after `flattenEnd`.

## Independent regression oracle

`tests/test_schedule_axis_pressure.py` does not reuse the production factorization/plateau antiderivatives as its main pressure oracle. It sends the unfactorized `final_angular(data,y,eta)^2` directly to SciPy adaptive improper quadrature on independently split intervals and compares that integral with the production result. Additional checks cover the official ideal-prefix lower pressure bound, evenness, derivative oddness/zero-at-origin, a finite-difference derivative cross-check, stable kernel evaluation, and fail-closed invalid inputs.

## Deliberate boundary

This closes the previous blocker **“materialize the actual `SchedulePressure.axisPressure` integral”** only at the level of an executable, reproducibly cross-checked numerical evaluator of the exact constructed schedule.

It does **not** provide a rigorous global quadrature error enclosure, a certified uniform positive minimum of `H^2` on the actual low-|Z| set, or the resulting paper-selected `sigma`. Sampled minima must not be promoted to that certificate. `Lambda/C`, the coefficient-space fixed-point fields, `NaturalProfileAssembly` on those real fields, and support/moment/matching/cone certificates remain unresolved. Therefore Stage 1 remains `formal-structure` and `paper_exact_velocity_available` remains false.
