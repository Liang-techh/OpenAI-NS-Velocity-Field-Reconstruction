# First-Picard axial `lin2(x1)` provenance

Status: **formal-structure only**.

This increment materializes the raw axial linear branch of the pinned
`naturalRemainder` at the genuine theorem-selected first Picard state:

```text
lin2(u) = j1((A*one - 4*A*(eta*uStar) + d*uStarEta) * u)
           + dot1(wStar, u)
           + param1(u, hStar)
```

The outer `inverseL` from the axial natural-remainder assembly is deliberately
not included here.

## Exact mixed-scale identity

The landed first Picard axial coefficient is represented as

```text
u1 = U0 + U1/Lambda + U2/Lambda^2 + a^2*P1/Lambda,
```

where `U0` is the reference coefficient, `U1` is the first ordinary
correction numerator, `U2` is the second ordinary correction numerator, and
`P1` is the normalized full pressure-derivative numerator.  The x1 state stores
`U1 = ordinary_base(x0)/2`, `U2 = -ordinary_slow(x0)/2`, and the pressure term
as the signed-log value `a^2*P1/Lambda`.

Because `lin2` is linear in its axial input, the returned coefficient is kept
as

```text
L(U0) + L(U1)/Lambda + L(U2)/Lambda^2 + a^2*L(P1)/Lambda,
```

with `L` equal to the three displayed `j1`, `dot1`, and `param1` terms.  The
negative slow sign is already present in `U2`; the `lin2` sum introduces no
additional minus sign.

The pressure source is read from
`wide_pressure.normalized_factor`, through
`ActualScheduleReferenceWideAxialRemainderState.pressure_normalized_factor`.
That value is the full pressure derivative jet divided by `a^2`, including the
amplitude derivatives.  The implementation divides it by 2 to obtain `P1`,
accumulates the normalized numerator in a local 96-digit Decimal context, and
rebuilds the final term as a `SignedLogCoefficientJet`.  It never differentiates
a zeroth-order normalized pressure value and never collapses a nonzero pressure
term into binary64.

## Pinned row action and domain

All fixed `AxisData` fields are radial degree zero and are rebuilt only from
`x1.reference`.  For `n = 0`, all three regular inverses return the exact zero
row.  For `n >= 1`, the source row is `n - 1`, the regular divisor is
`radialDivisor(1,n-1) = n^2`, `dot1` carries the Euler factor `n - 1`, and
`param1` reads the source at parameter order `k + 1`.  The executable domain is
`n,m >= 0` and `eta ∈ [-11/10, 11/10]`; invalid indices, eta values, state
types, epsilon values, or theorem scales fail closed.

The state accepts only a genuine `ActualScheduleWideFirstPicardState` with
`picard_x1_materialized = true`, a positive theorem-selected `Lambda`, and the
matching actual pressure-amplitude scale.  It exposes
`lin2_x1_materialized = true` and `axial_lin2_materialized = true`, while
`natural_remainder_x1_materialized`, `fixed_point_materialized`,
`fixed_point_convergence_certified`, `global_axis_norm_certified`, and
`paper_exact` remain false.

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
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_lin2.py`.
The upstream expression is mirrored by the existing raw remainder at
`src/openai_ns_reconstruction/axis_coefficient_natural_remainder.py:151-157,195-201`;
the actual x0 axial split is established at
`src/openai_ns_reconstruction/axis_coefficient_wide_axial_remainder.py:297-334`,
and the x1 split is established at
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard.py:149-187,247-282`.

## Independent verification boundary

The local test file passed 8 tests, including direct production evaluation on
diagnostic polynomials with independently expanded rational reference values.
Measured scope and limitations are recorded in
`reports/stage1_first_picard_local_integration.md`.

An independent test should evaluate the three fixed operator formulas on the
reference, first ordinary numerator, and second ordinary numerator separately,
then compare the returned Decimal numerators.  It should independently repeat
the radial row action with source row `n-1`, divisor `n^2`, Euler factor `n-1`,
and the `k+1` parameter derivative.  For pressure, the test should reconstruct
the normalized full derivative source from `pressure_normalized_factor`, divide
by 2, and compare the signed-log sign, common amplitude log, and
`log(abs(numerator)) - log(Lambda)` factor.  Row zero, invalid domains,
mismatched/non-x1 states, and nonzero pressure that cannot fit binary64 should
remain fail-closed.

The companion `slow2` module now supplies its complete aggregate. This adapter
does not perform their outer recombination or materialize `lin1`, `quad1`, `slow1`,
the x1 source/pressure/resolvent path, `naturalRemainder(x1)`, a later Picard
iterate, convergence, final `phi/u`, derived average/pressure fields,
`NaturalProfileAssembly`, global weighted `AxisSpace` bounds, or paper-exact
velocity data.
