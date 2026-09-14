# Axis coefficient wide first-Picard `lin2(x1)` provenance

This increment is a narrow Stage-1 backend artifact stacked on PR #352. It
materializes only the pinned axial-linear constituent of `naturalRemainder` at
the already-materialized genuine actual-schedule first Picard state.

Pinned formal source:

- repository: `openai/NavierStokesAndEuler`
- revision: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisContraction.lean`:
  `axialLinearCoefficient` and the `lin2` term inside `naturalRemainder`
- `NavierStokes/AxisOperators.lean` / `AxisWeightEstimates.lean`:
  coefficient product, zero-datum `regularInverse` with `r=1`,
  `inverseDotProduct`, `inverseParamProduct`, and
  `radialDivisor r n = (n+1)(n+r)`

The implemented identity is exactly

`lin2(x1) = j1(product(axialLinearCoefficient(d),u1))
            + dot1(wStar,u1) + param1(u1,hStar)`

with

`axialLinearCoefficient(d) =
 A*one - 4*A*product(eta,uStar) + product(d.d,uStarEta)`.

For output radial row `n>0`, all three branches use the pinned `r=1`
zero-datum inverse divisor `radialDivisor(1,n-1)=n^2`. `dot1` applies the
radial Euler factor to its second argument, and `param1` takes one eta
derivative of its first argument before the product.

The implementation consumes `ActualScheduleWideFirstPicardState` directly.
It preserves the ordinary `u1` families `Lambda^0`, `Lambda^-1`,
`Lambda^-2` and separately carries the genuine pressure-derived
`a(eta)^2/Lambda` family via the actual schedule pressure normalized factor
and amplitude log. No caller-supplied coefficient table, amplitude, fitted
field, sampled identity, or absent-row default is accepted.

Truth boundary: this proves/materializes only this finite constituent at x1.
`slow1(x1)`, the complete `naturalRemainder(x1)`, x2, a converged fixed point,
global `AxisCoefficientSpace` norm certification, `NaturalProfileAssembly`,
and paper-exact velocity remain unavailable.
