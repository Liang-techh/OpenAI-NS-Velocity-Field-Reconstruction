# Stage 1 validated phase local integration

Updated 2026-09-13. This report records the finite exact rational phase seam and its actual schedule
amplitude bridge. It is a formal structure component, not a paper exact amplitude or reconstruction.

## Pinned relation

At pinned commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
`NaturalAxisData.lean:24-35` gives

```text
L(h,x) = 1 - 2 h x^2
H(h,j,x) = (1/2-h)x + (1-x^2)(4x+j)
```

`NaturalAxisCoefficients.lean:47-49, 413-414, 468-485` gives

```text
g(x) = -L(h,x) H(h,j,x)/(H(h,j,x)^2+sigma^2)
phase(eta) = integral_0^eta g(x) dx
a(eta) = exp(Lambda phase(eta))/C
```

The negative eta orientation is handled by integrating `[eta,0]` and negating the result. The
coefficient window is `[-11/10,11/10]` (`NaturalAxisCoefficients.lean:23-27`).

## Implemented seam

`src/openai_ns_reconstruction/axis_phase_integral.py` expands `P=-L H` and
`Q=H^2+sigma^2` around each rational cell center. It accepts a cell only when

```text
theta = sum(k>=1) abs(Q_k) (2r)^k / Q_0 < 1/2
```

and then uses the exact rational Taylor recurrence and Cauchy tail. Cell budgets sum to the requested
phase tolerance; midpoint subdivision and Taylor order growth are bounded by finite caps and fail
closed. Refinement at fixed order is not a termination argument, so the implementation raises order
as well as subdividing. `eta=0` is an exact zero phase case.

`ActualScheduleAmplitudeLogState.log_amplitude_enclosure` in
`src/openai_ns_reconstruction/axis_coefficient_amplitude.py:225-272` converts the selected actual
`h,j,sigma,eta` inputs to the rational path and requests
`absolute_log_tolerance/Fraction(Lambda)`. Since `Lambda` is a finite positive `Decimal`, its Fraction
conversion is exact. The returned interval is transported through
`log(a)=Lambda*phase-log(C)` with directed Decimal presentation.

The executable scale selector stores `C` symbolically as `exp(C_exponent)`
(`natural_scale_selection_wide.py:93-110,169-191`), where `C_exponent` is the upward rounded selected
finite exponent. The bridge therefore treats that stored exponent as the selected `log(C)`; it does
not claim an enclosure for an unspecified lower value of `log(C)`.

## Actual fixture and resource evidence

The shared fixture uses `P=2`, `m=1`, `lam=0.05`, `wait=30`, `h=0.01`, `j=0.05`. Its actual margin path
selects `sigma = 1.518624346196426e-9`; the formal root lies in `(-0.0125,-0.01)` and the numerical
root is about `-0.011135706745859649`. A direct wide chain inspection gives
`Lambda = 4.774138249520899e+775` and `C_exponent = 4.722378862791228e+784`.

The focused command
`python -m pytest -q tests/test_axis_phase_integral.py -W error` returned `5 passed in 0.17 s`.
The actual amplitude bridge command
`python -m pytest -q tests/test_axis_amplitude_phase_enclosure.py -W error` returned `3 passed in
0.17 s`; the existing amplitude compatibility command returned `6 passed in 0.16 s`.

An actual root crossing probe at `eta=-1/50` with coarse `absolute_log_tolerance=10^777` returned in
`1.559 s` through the amplitude bridge and `1.473 s` through the direct phase integrator. It accepted
`225` cells, used maximum Taylor order `20`, and reported phase error within the phase budget and log
interval width within `2*10^777`. This is coarse resource evidence only. A strict high order estimate
of roughly `2,600` Taylor orders for a much smaller log tolerance is not a benchmark.

The same actual selected binary64 inputs were then run at direct phase tolerance `10^-12`. The exact
integrator returned in `10.664 s`, with `225` accepted cells, maximum order `64`, estimate
`0.05073499503328993`, and error bound `7.357907e-13` (the bound is below `10^-12`). The legacy
4001 point `real_phase` call returned in `0.002 s` with value `1.521550417558555`; its
`Fraction.from_float` value lies outside the validated interval by approximately `+1.470815`.
This is one actual selected data fixture, not a general benchmark, and the exact rational layer and
legacy floating path should not be described as identical arithmetic.

## Boundary

The default 4001 point `real_phase` trapezoid remains the point estimate used by `log_amplitude`, and
the measurement above makes replacing or reworking that default phase path the next blocker for
crossing-root evaluations. The validated seam covers only the selected rational kernel and affine log
amplitude transport. It does not certify SchedulePressure quadrature, theorem parameter selection,
coefficient roundoff, global weighted `AxisSpace` membership, or the full reconstruction. Finite cell
and order caps fail closed; no global approximation or convergence claim follows from these checks.
