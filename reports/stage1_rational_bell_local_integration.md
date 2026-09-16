# Stage 1 rational Bell local integration

Status: **implemented; affected-source and integration acceptance passed**.

`RationalAxisCoefficientData` now supplies
`amplitude_power_bell_fraction(power, order, eta, Lambda)`. It interprets the
selected `Lambda` Decimal as an exact rational and forms

```text
q_k = power * Lambda * jet_fraction("normalizedGradient", 0, k, eta)
B_0 = 1
B_(n+1) = sum_k binom(n,k) q_k B_(n-k)
```

over exact `Fraction` values. The result is
`(d_eta^order a^power) / a^power`; no amplitude, phase, or `C` is evaluated.
`amplitude_power_bell_decimal` uses explicit nearest Decimal96 rounding, and
`amplitude_power_bell_enclosure` returns immutable directed Decimal96 bounds.

The input checks require nonnegative integer `power` and `order`, positive
finite Decimal `Lambda`, and eta in the exact window `[-11/10, 11/10]`.
The selected parameters remain exact binary rationals; this is not pressure,
phase, theorem-selection, global `AxisSpace`, or accumulated-error certification.

The affected-source acceptance command was:

```text
python -m pytest -q tests\test_axis_coefficient_wide_natural_source.py tests\test_axis_coefficient_rational_bell.py -W error
```

It returned `12 passed in 0.30 s`, after replacing old float-gradient
expectations with independent Fraction expectations for the first two Bell
values. The broader formal/Picard/profile/wide integration command

```text
python -m pytest -q tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_picard_family.py tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error
```

returned `13 passed in 24.37 s`; phase and total-product construction and
accumulated coefficient-error control remain outside this seam.

The bounded source check passed:

```text
python -B -m py_compile src/openai_ns_reconstruction/axis_coefficient_rational_data.py
```

No additional pytest or broader validation was run after these reported checks;
the remaining cleanup is compile-only.
