# First-Picard angular nonlinear provenance

Status: **formal-structure only**.

Artifacts:

- `src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_angular_nonlinear.py`
- `references/AXIS_COEFFICIENT_WIDE_FIRST_PICARD_ANGULAR_NONLINEAR_PROVENANCE.md`
- `references/provenance_manifest_addendum_axis_coefficient_wide_first_picard_angular_nonlinear.json`

This increment evaluates both pinned angular nonlinear expressions at the
same actual theorem-selected first Picard state `x1`:

```text
quad1 = j2((d * normalizedGradient) * u1 * phi1)

slow1 = j2((averageCoefficient * average(u1)
            + angularSlowCoefficient * u1) * phi1)
        + param2(average(u1), d * phi1)
        + dot2(averageCoefficient * average(u1), phi1)
        + d * mixed2(average(u1), phi1)
        - param2(phi1, d * u1)
```

The factory is `wide_first_picard_angular_nonlinear_state(x1)`.  Its state
exposes `quad1_jet(n, m, eta)` and `slow1_jet(n, m, eta)`, each returning a
frozen `MixedScaleFirstPicardAngularNonlinearCoefficientJet` with:

- `ordinary`: five Decimal numerators for `Lambda^0` through `Lambda^-4`;
- `pressure_linear`: three normalized Decimal numerators for
  `a^2 Lambda^-1` through `a^2 Lambda^-3`;
- the actual `Lambda` and the actual amplitude logarithm.

All branches are rebuilt from the actual typed `x1` and its
`ActualScheduleAxisCoefficientData`.  Radial products use the finite
`i + j = n` convolution and parameter products use the binomial eta
Leibniz sum.  `j2`, `param2`, `dot2`, and `mixed2` use the exact divisor
`n * (n + 1)`.  The pressure channel is supplied by the actual axial
`pressure_normalized_factor / 2`; pressure-pressure products and any nonzero
term outside the declared channel widths fail closed, so no `a^4` term is
silently discarded.  The final pinned minus sign on the last `param2` term is
applied exactly once.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Modules:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/AxisWeightEstimates.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`

The exact upstream coefficient definitions are in
`AxisContraction.lean:327-334`, and the `naturalRemainder` source assignments,
including `quad1` and `slow1`, are in `AxisContraction.lean:339-365`.  The
operator laws used here are
`AxisOperators.lean:86-95` (`product`), `AxisOperators.lean:501-504`
(`regularInverse`), `AxisOperators.lean:512-518` (differential-family
operators), and `AxisOperators.lean:681-717` (coefficient bindings).  The
radial divisor is pinned at `AxisWeightEstimates.lean:388`; the coefficient
window is pinned at `NaturalAxisCoefficients.lean:27`.

The local executable source map is
`src/openai_ns_reconstruction/axis_coefficient_natural_remainder.py:140-193`
for the fixed angular coefficients and both expressions, and
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard.py:247-282`
for the mixed-scale `x1` pair.  The axial pressure normalization comes from
`src/openai_ns_reconstruction/axis_coefficient_wide_axial_remainder.py`.

## Numerical and truth boundary

The channel arithmetic uses a local 96-digit Decimal context.  It is a
numerical representation of the pinned coefficient formulas, not exact real
arithmetic, an interval enclosure, or a theorem certificate.  Rounded
cancellation can lose digits; the implementation keeps the formal-structure
truth status and does not promote the result.  Signed-log pressure views keep
the common amplitude log separate from the finite coefficient and Lambda
factors, and no binary64 conversion of the pressure channel is performed.

This layer does not apply the outer `inverseL` or natural resolvent.  It does
not materialize `naturalRemainder(x1)`, a later Picard iterate, a fixed point,
convergence, a global weighted `AxisSpace` certification, a derived final
profile, or a paper-exact velocity.  The phase and amplitude logarithm retain
the landed actual-schedule numerical quadrature boundary.

The state therefore reports `paper_exact=False`,
`global_axis_norm_certified=False`, `natural_remainder_x1_materialized=False`,
`fixed_point_materialized=False`, and
`fixed_point_convergence_certified=False`.
