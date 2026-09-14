# Genuine first-Picard `quad1(x1)` provenance

## Scope

This increment repairs and materializes exactly the pinned angular quadratic branch

`quad1(x1) = j2(product(product(angularQuadraticCoefficient(d), u1), phi1))`

with

`angularQuadraticCoefficient(d) = product(d.d, d.normalizedGradient)`.

The earlier branch formula `-j2(product(primitive(phi1),phi1))` was not the
pinned `AxisContraction.naturalRemainder` definition and is removed rather than
carried forward as a compatible approximation.

The landed genuine `x1` is used in full. `phi1` contains ordinary
`Lambda^0, Lambda^-1, Lambda^-2` families. `u1` contains those same ordinary
families plus the genuine pressure-derived `a^2/(2*Lambda)` contribution. The
ordinary `quad1` result is therefore retained through `Lambda^-4`, and the
pressure-derived part is retained as `a^2` times `Lambda^-1, Lambda^-2,
Lambda^-3` normalized factors. No binary64 collapse of the theorem-scale
pressure amplitude is performed.

This does **not** claim that `naturalRemainder(x1)` is complete, does not produce
`x2`, and does not certify global `AxisCoefficientSpace` membership,
contraction closure, a fixed point, `NaturalProfileAssembly`, or paper-exact
velocity.

## Formal source

Pinned source revision:

- repository: `openai/NavierStokesAndEuler`
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisContraction.lean`
  - `NavierStokes.AxisContraction.angularQuadraticCoefficient`
  - `NavierStokes.AxisContraction.naturalRemainder`
  - literal branch `quad1 := O.j2 (O.product (O.product (angularQuadraticCoefficient O d) u) phi)`
- pinned coefficient operators:
  - `product` uses finite radial convolution and the eta Leibniz rule;
  - `j2` is `regularInverse(_,2)`, hence output row zero is literal zero and
    output row `n>0` reads input row `n-1` divided by `n(n+1)`.

There is no minus sign and no `primitive(phi)` in the pinned `quad1` branch.

## Actual-data dependency

The production factory accepts only the landed
`ActualScheduleWideFirstPicardState`. From that same source-bound object it
consumes:

- actual `AxisData.d` and `AxisData.normalizedGradient`;
- the landed `coefficientOperators.product` implementation;
- genuine mixed-scale `u1` and `phi1` jets;
- the actual SchedulePressure wide pressure normalized factors and amplitude;
- the theorem-selected `Lambda` and common coefficient epsilon.

The wide axial pressure factor is not dropped. If `p[n,m]` is the landed
remainder pressure derivative after factoring out `a^2`, the x1 pressure piece
is represented exactly as `a^2 * p[n,m] / (2*Lambda)` before the remaining
ordinary products and `j2` map are applied.

No caller-supplied coefficient table, cutoff, Lambda, amplitude, surrogate
field, fitted sample, or absent-row default is accepted.

## Regression obligations

`tests/test_axis_coefficient_wide_first_picard_quad1.py` checks:

1. binding to genuine x1, actual AxisData, and the genuine wide-pressure chain;
2. literal `j2` row-zero behavior for both ordinary and pressure scale families;
3. the first radial row identity for `q*u1*phi1/2`, including the pressure part;
4. a general `(n,m)` coefficient against an independently replayed triple radial
   convolution, three-factor eta Leibniz sum, pressure factor, and `j2` divisor;
5. signed-log retention of the `a^2 Lambda^-1..-3` pressure-derived families;
6. rejection of surrogate states and invalid indices/window coordinates; and
7. machine-readable provenance with `naturalRemainder(x1)`, `x2`, fixed-point,
   and paper-exact claims still false.

## Remaining boundary

The genuine x1 source, x1 pressure, complete four-constituent `slow2(x1)`, and
now the corrected pinned `quad1(x1)` are available on compatible actual-schedule
data. The shortest remaining closure path is to materialize the pinned
`lin1(x1)`, `slow1(x1)`, and `lin2(x1)` mixed-scale branches and then assemble
the complete angular/axial `naturalRemainder(x1)`. Only that complete typed
pair may be passed through the outer Picard map to construct genuine `x2`.
Iteration, contraction/tail closure, global coefficient-space membership,
fixed-point `phi/u`, derived average/pressure fields, `NaturalProfileAssembly`,
and paper-exact velocity remain unresolved.
