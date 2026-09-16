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

The accumulation baseline keeps the public exact-Fraction integral mode for
compatibility while routing the default amplitude through an outward dyadic
cell accumulator. It uses exact local `[I-E,I+E]` intervals, a dyadic quantum
`d <= b/4` for cell budget `b`, and floor/ceiling endpoint projection. With
`E <= b/2`, each projected cell radius is at most `E+d <= 3b/4`, so the summed
radius is at most `3/4` of the requested tolerance. Dyadic endpoints sum without
a general-rational LCM. This is an implementation contract; it does not weaken
the requested tolerance or change the selected phase inputs.

## Evidence

The historical pre-default focused commands returned `5 passed in 0.17 s` for
`tests/test_axis_phase_integral.py`, `3 passed in 0.17 s` for
`tests/test_axis_amplitude_phase_enclosure.py`, and `6 passed in 0.16 s` for
the prior amplitude compatibility target. These `14` relevant tests predate the
default-path replacement and remain historical evidence.

The earlier checkpoint command combining
`tests/test_axis_amplitude_validated_default.py`,
`tests/test_axis_coefficient_amplitude.py`, and
`tests/test_axis_phase_log_enclosure.py` was interrupted after approximately
`9` minutes of CPU-heavy exact phase work. It produced no pytest summary or
failure traceback, so no pass or failure result is claimed.

A bounded diagnostic at `eta=-1/50` completed the rational phase in `10.850 s`
with `225` cells and maximum order `64`; affine log transport took `0.558 s`.
The Decimal midpoint path measured `11.307 s`, with stack evidence in
Decimal/Fraction conversion. The current default validated binary64-input
diagnostic began at `23.019 s` and hit the `45 s` bound before completing; no
precise runtime is claimed.

The completed slow2 validation recorded
`python -m pytest -q tests/test_axis_phase_integral.py -W error` as `7 passed in
0.20 s`. The combined command
`python -m pytest -q tests/test_axis_amplitude_validated_default.py
tests/test_axis_coefficient_amplitude.py tests/test_axis_phase_log_enclosure.py
-W error` returned `13 passed in 43.71 s` under the hard `120 s` cap. This
resolves the default-regression performance blocker for that bounded command.
The earlier approximately `9` minute interruption remains historical and had no
pytest summary; no full-suite or paper-exact claim follows.

The amplitude log-scale foundation command
`python -m pytest -q tests/test_axis_amplitude_log_scale.py -W error` reports
`5 passed in 0.15 s`. Its focused integration command
`python -m pytest -q tests/test_axis_amplitude_log_scale_integration.py -W error`
reports `6 passed in 0.20 s`. These checks cover actual eta-zero amplitude,
source, pressure, and axial rows, a synthetic q2 helper, and scale-only copies;
they do not establish full coefficient arithmetic.
The mixed-route command
`python -m pytest -q tests/test_axis_amplitude_log_scale_mixed.py -W error`
reports `9 passed in 0.22 s`; the aggregate command
`python -m pytest -q tests/test_axis_amplitude_log_scale_aggregate.py -W error`
reports `4 passed in 0.18 s`. These four focused commands total `24` distinct
tests. The four slow2 branches and aggregate path anchor the actual source and
check cross-branch compatibility before summing, including the zero branch.
Parameter and average-dot mixed q2/q4 routes are implemented, compiled, and
root-reviewed. The metadata chain now reaches all assigned nonlinear branches,
the complete coefficientwise remainder/second-Picard path, the formal triangular
solver, and the wide-profile surfaces. Full baseline integration remains pending.
The full-suite run immediately preceding this metadata work returned
`1936 passed, 1 failed` in `1723.45 s`; the sole fixture-argument issue was fixed
and its affected test subsequently passed in `0.20 s`, but the full suite has not
been rerun after that fix.

The existing selected-fixture measurement at `eta=-1/50` used phase tolerance
`10^-12`, accepted `225` cells, reached Taylor order `64`, and returned phase
estimate `0.05073499503328993` with error bound `7.357907e-13` in `10.664 s`.
The current default validated binary64-input diagnostic hit its `45 s` cap before
completing, so no precise runtime is claimed. These are selected-fixture
observations, not a general performance benchmark.

## Boundary and next target

The interval certifies the chosen rational phase kernel and the affine transport
through the selected finite `Lambda` and symbolic `C` exponent. It does not
certify SchedulePressure quadrature, theorem parameter selection, `C` selection,
coefficient roundoff, global weighted `AxisSpace` membership, or the full
reconstruction. The already-landed coefficientwise `naturalRemainder(x0)`,
`x1`, `naturalRemainder(x1)`, `x2`, formal Picard family, profile prefix, and
conditional profile-budget interfaces retain their documented formal/conditional
status. The next Stage 1 boundary is exact Bell-factor/logarithm enclosures,
general coefficient arithmetic, SchedulePressure and parameter proofs, and
certifying compatible global `AxisSpace` norms; no global certification follows
from this default phase correction.

## Current checkpoint — 2026-09-14

The metadata evidence totals `35` distinct tests: prior `24`, plus `5` nonlinear
tests in `0.80 s`, `3` remainder/second-Picard tests in `0.38 s`, and `3`
formal/wide-profile tests in `0.22 s`. A separate batch of the existing default
test file returned `4 passed in 20.19 s` and included one new nonzero-eta
(`eta=-0.02`) metadata test, bringing the new metadata total to `36` distinct
tests. The demo exited `0` with all checks. The direct strict-audit module exited
`2` as expected and this was confirmed in its log; an earlier console-wrapper
exit `1` is historical. All relevant source compilation and review checks passed.

The local checkpoint is being saved without an external push; current full-suite
acceptance remains pending.

