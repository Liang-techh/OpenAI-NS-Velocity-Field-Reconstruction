# First-Picard angular `lin1(x1)` provenance

Status: **formal-structure only**.

This increment materializes the raw angular linear branch of the pinned
`naturalRemainder` at the genuine theorem-selected first Picard state:

```text
lin1(phi) = j2((wStar + h*one - 2*h*(eta*uStar)) * phi)
             + dot2(wStar, phi)
             + param2(phi, hStar)
```

The returned object reuses the existing
`MixedScaleFirstPicardAngularCoefficientJet` split representation.  The outer
`inverseL` multiplication and angular `naturalResolvent` are deliberately not
included here.

## Exact mixed-scale identity

The first Picard angular coefficient is represented as

```text
phi1 = Phi0 + Phi1/Lambda + Phi2/Lambda^2.
```

Because `lin1` is linear in its angular input, the returned channels are

```text
L(Phi0) + L(Phi1)/Lambda + L(Phi2)/Lambda^2,
```

where `L` is the displayed `j2`, `dot2`, and `param2` sum.  No signed-log
pressure channel is present in this angular branch.

## Pinned row action and domain

All fixed `AxisData` fields are radial degree zero and are rebuilt only from
`x1.reference`.  For `n = 0`, all three zero-datum regular inverses return the
exact zero row.  For `n >= 1`, the source row is `n - 1`, the regular divisor is
`radialDivisor(2,n-1) = n*(n+1)`, `dot2` carries the Euler factor `n - 1`, and
`param2` reads the angular source at parameter order `k + 1`.  The executable
domain is `n,m >= 0` and `eta ∈ [-11/10, 11/10]`; invalid indices, eta values,
state types, epsilon values, or theorem scales fail closed.

The state accepts only a genuine `ActualScheduleWideFirstPicardState` with
`picard_x1_materialized = true`, a positive theorem-selected `Lambda`, and
matching actual coefficient epsilon.  It exposes
`lin1_x1_materialized = true`, `angular_lin1_materialized = true`, and
`lin1_complete = true`, while `natural_remainder_x1_materialized`,
`fixed_point_materialized`, `fixed_point_convergence_certified`,
`global_axis_norm_certified`, and `paper_exact` remain false.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Files:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/AxisWeightEstimates.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`
- Symbols:
  - `NavierStokes.AxisContraction.naturalRemainder`
  - `NavierStokes.AxisOperators.productFamily`
  - `NavierStokes.AxisOperators.regularInverse`
  - `NavierStokes.AxisOperators.inverseDotProduct`
  - `NavierStokes.AxisOperators.inverseParamProduct`
  - `NavierStokes.NaturalAxisCoefficients.CoefficientFamily.axisData`

The Python source is
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_lin1.py`.
The pinned raw expression is mirrored by
`src/openai_ns_reconstruction/axis_coefficient_natural_remainder.py:129-170`;
the first Picard angular split is established by
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard.py:105-145,247-282`.

## Independent verification boundary

An independent test should evaluate the three fixed operator formulas on the
reference, first ordinary numerator, and second ordinary numerator separately,
then compare the returned Decimal channels.  It should independently repeat
the radial row action with source row `n - 1`, divisor `n*(n+1)`, Euler factor
`n - 1`, and the `k + 1` parameter derivative.  Row zero, invalid domains,
mismatched/non-x1 states, and all downstream truth flags should remain
fail-closed.

This adapter does not materialize the angular `quad1` or `slow1` terms, the
axial `lin2` or `slow2` terms, the x1 source/pressure/resolvent path, complete
`naturalRemainder(x1)`, a later Picard iterate, convergence, final `phi/u`,
derived average/pressure fields, `NaturalProfileAssembly`, global weighted
`AxisSpace` bounds, or paper-exact velocity data.
