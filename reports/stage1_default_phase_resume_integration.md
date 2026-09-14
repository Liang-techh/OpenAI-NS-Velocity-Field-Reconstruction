# Stage 1 default phase resume integration

Updated 2026-09-14. This report records the default phase-path correction after
the validated exact-rational phase seam. It is a conditional scalar-input
component, not a paper-exact amplitude or reconstruction.

## Default path

`ActualScheduleAmplitudeLogState` now uses the validated Fraction phase integral
for its default `log_amplitude` value. The default phase tolerance is
`1/10^12`; exact selected `h,j,sigma,eta` inputs and finite integrator settings
are cached with capacity `16`. The default value is the nearest Decimal96
midpoint of the returned log-amplitude interval.

`default_log_amplitude_enclosure` exposes the phase interval after the selected
positive `Lambda` amplification and symbolic `C` exponent transport. The
explicit `log_amplitude_enclosure` method accepts a rational log tolerance and
passes the exact quotient `absolute_log_tolerance/Fraction(Lambda)` to the
phase integrator. Finite cell, depth, and Taylor-order caps fail closed.

`phase_samples` is retained only for
`legacy_log_amplitude_diagnostic`, which names the former 4001-point trapezoid.
That path is diagnostic and no longer supplies the default amplitude value.

## Evidence

The existing focused commands returned `5 passed in 0.17 s` for
`tests/test_axis_phase_integral.py`, `3 passed in 0.17 s` for
`tests/test_axis_amplitude_phase_enclosure.py`, and `6 passed in 0.16 s` for
the prior amplitude compatibility target. These `14` relevant earlier tests
are historical evidence from before the default-path replacement and do not
accept the new default.

The checkpoint command
`python -m pytest -q tests/test_axis_amplitude_validated_default.py
tests/test_axis_coefficient_amplitude.py tests/test_axis_phase_log_enclosure.py
-W error` was interrupted after approximately `9` minutes of CPU-heavy exact
phase work. It produced no pytest summary or failure traceback, so no pass or
failure result is claimed. The new default remains implemented but unaccepted
until its performance and regression evidence are resolved; no additional run
was made in this batch.

The existing selected-fixture measurement at `eta=-1/50` used phase tolerance
`10^-12`, accepted `225` cells, reached Taylor order `64`, and returned phase
estimate `0.05073499503328993` with error bound `7.357907e-13` in `10.664 s`.
The legacy 4001-point diagnostic returned `1.521550417558555` in `0.002 s` and
lay outside the validated interval by approximately `+1.470815`. This is one
selected fixture, not a general performance benchmark; the exact rational and
legacy floating paths are not identical arithmetic.

## Boundary and next target

The interval certifies the chosen rational phase kernel and the affine transport
through the selected finite `Lambda` and symbolic `C` exponent. It does not
certify SchedulePressure quadrature, theorem parameter selection, `C` selection,
coefficient roundoff, global weighted `AxisSpace` membership, or the full
reconstruction. The already-landed coefficientwise `naturalRemainder(x0)`,
`x1`, `naturalRemainder(x1)`, `x2`, formal Picard family, profile prefix, and
conditional profile-budget interfaces retain their documented formal/conditional
status. The next Stage 1 boundary is propagating scalar-input and coefficient
error bounds through the genuine wide natural-remainder/profile chain and then
certifying compatible global `AxisSpace` norms; no global certification follows
from this default phase correction.

Merge resolution: all five add/add conflicts were resolved and staged, preserving the local complete aggregate and the new-main parameter branch. Source compilation and both JSON parses passed. The focused slow2 batch reported three completed tests over 318.78 seconds before interruption in the phase path; this is partial evidence, not complete acceptance. Latest main integrated locally: 065b67e. No remote push is claimed.

