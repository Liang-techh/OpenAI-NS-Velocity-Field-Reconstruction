# Stage 1 rational-input local integration

Status: **implemented; measured baseline passed**.

`axis_coefficient_rational_data.py` adds `RationalAxisCoefficientData` anchored
only by `ActualScheduleAxisCoefficientData`. It interprets selected `h`, `j`,
and `sigma` values as exact binary rationals, derives `A = 1/2 + h` and
`D = 1/2 - h`, and evaluates the nine fixed scalar fields plus `chi` with
normalized Taylor algebra over `Fraction`.

The provider exposes `jet_fraction`, `angular_reference_jet_fraction`, nearest
Decimal96 `jet_decimal` and `angular_reference_jet_decimal`, and immutable
directed `RationalJetEnclosure` values, including an angular enclosure wrapper.
The exact input window is `[-11/10, 11/10]`; binary64 literal `1.1` is rejected
because it is slightly outside exact `11/10`, while the physical `|eta| <= 1`
chart is unchanged. `zStar` is explicitly unsupported and is never replaced by
zero.

This seam does not provide pressure/Bell-phase inputs, a theorem parameter
selection proof, global `AxisSpace` membership, or accumulated coefficient-error
control. The measured baseline command

```text
python -m pytest -q tests\test_axis_coefficient_rational_data.py tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_picard_family.py tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error
```

returned `16 passed in 25.22 s`. The remaining cleanup only removes a discarded
reference-half evaluation and preserves the value formulas.

The bounded source check passed:

```text
python -B -m py_compile src/openai_ns_reconstruction/axis_coefficient_rational_data.py
```

No additional pytest or broader validation was run after this measured baseline;
the remaining cleanup is compile-only.
