# Complete first-Picard `slow2` aggregate provenance

Status: **formal-structure only**.

Artifacts:

- `src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_slow2.py`
- `tests/test_axis_coefficient_wide_first_picard_slow2.py`

This increment assembles the four pinned terms in
`AxisContraction.naturalRemainder.slow2` at the genuine theorem-selected first
Picard state `x1`:

```text
j1((2*A*eta) * (u1*u1))
  + dot1((2*D*eta) * average(u1), u1)
  + d * mixed1(average(u1), u1)
  - param1(u1, d*u1)
```

The four branch states are evaluated from the same actual `x1`, then their
ordinary and normalized pressure numerators are summed coefficientwise.  The
ordinary part stays split over `Lambda^0` through `Lambda^-4`; the pressure
linear part stays normalized over `a^2 Lambda^-1` through `a^2 Lambda^-3`; and
the pressure-square part stays normalized over `a^4 Lambda^-2`.  The parameter
branch already includes its displayed minus sign, so aggregation adds that
branch as supplied.

The aggregate uses ordinary Decimal arithmetic in a local 96-digit context.
Those sums are a numerical representation of the pinned coefficient formulas,
not exact real arithmetic, an interval enclosure, or an independent theorem
certificate.  Rounded cancellation in a summed numerator can lose digits; the
implementation does not promote that result to a stronger truth status.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main files:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/AxisWeightEstimates.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`
- Relevant symbols:
  - `NavierStokes.AxisContraction.naturalRemainder`
  - `NavierStokes.AxisOperators.regularInverse`
  - `NavierStokes.AxisOperators.inverseDotProduct`
  - `NavierStokes.AxisOperators.inverseMixed`
  - `NavierStokes.AxisOperators.inverseParamProduct`
  - `NavierStokes.AxisOperators.average`
  - `NavierStokes.NaturalAxisCoefficients.CoefficientFamily.axisData`

The exact upstream source locations used for this boundary are
`AxisContraction.lean:327-334` for the coefficient definitions and
`AxisContraction.lean:339-365` for `naturalRemainder`; the operator bindings are
`AxisOperators.lean:86-95` (`product`), `AxisOperators.lean:501-504`
(`regularInverse`), `AxisOperators.lean:512-518` (`differentialFamily`), and
`AxisOperators.lean:681-717` (the coefficient-operator record bindings).
The radial divisor is pinned at `AxisWeightEstimates.lean:388`, and the
coefficient window is pinned at `NaturalAxisCoefficients.lean:27`.

The aggregate uses the branch identities already pinned and individually
tested in the four companion modules.  It does not introduce a radial or
parameter truncation and does not convert a theorem-scale pressure term to
binary64.

## Repository boundary

The aggregate is based on current repository main commit
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.  The parameter branch is supplied
by the reviewed open PR 279 prerequisite, whose head is
`11811fb66101bf31b5be892062f6291f3091a410` and whose four files are:

- `src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_slow2_param.py`
- `tests/test_axis_coefficient_wide_first_picard_slow2_param.py`
- `references/AXIS_COEFFICIENT_WIDE_FIRST_PICARD_SLOW2_PARAM_PROVENANCE.md`
- `references/provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_param.json`

## Truth boundary

`ActualScheduleWideFirstPicardSlow2State` exposes
`slow2_complete=True` only after constructing and validating all four branch
states from the same `ActualScheduleWideFirstPicardState`, including shared
`epsilon`, `Lambda`, and per-jet amplitude logarithm.  It keeps
`natural_remainder_x1_materialized=False`, `fixed_point_materialized=False`,
`fixed_point_convergence_certified=False`, `paper_exact=False`, and the global
weighted `AxisSpace` certification false.

The final aggregate test file passed 9 tests, including nonzero pressure,
cancellation, and amplitude-identity checks on the production aggregation path.
See `reports/stage1_first_picard_local_integration.md` for measured limitations.

The companion module now materializes raw axial `lin2(x1)`. This aggregate does
not perform their outer recombination or materialize angular `lin1`, `quad1`, or
`slow1`, the `x1` source/pressure/resolvent path,
`naturalRemainder(x1)`, a later Picard iterate, convergence to the fixed point,
derived final average/pressure fields, `NaturalProfileAssembly`, or the
downstream support, moment, matching, cone, and paper-exact velocity claims.
