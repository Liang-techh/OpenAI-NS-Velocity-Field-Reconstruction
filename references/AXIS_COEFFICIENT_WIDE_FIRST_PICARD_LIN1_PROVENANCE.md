# Genuine first-Picard `lin1(x1)` provenance

## Scope

This increment materializes exactly the pinned angular-linear constituent

`lin1(x1) = j2(product(angularLinearCoefficient(d), phi1)) + dot2(wStar, phi1) + param2(phi1, hStar)`

with

`angularLinearCoefficient(d) = wStar + h*one - 2*h*product(eta,uStar)`.

It consumes the already-materialized genuine actual-schedule first Picard state and preserves the three ordinary `phi1` families `Lambda^0`, `Lambda^-1`, and `Lambda^-2` independently. No pressure-derived scale belongs to this branch because pinned `lin1` depends only on `phi1` and fixed AxisData.

This does **not** claim that `naturalRemainder(x1)` is complete, does not produce `x2`, and does not certify global `AxisCoefficientSpace` membership, contraction closure, a fixed point, `NaturalProfileAssembly`, or paper-exact velocity.

## Formal source

Pinned source revision:

- repository: `openai/NavierStokesAndEuler`
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisContraction.lean`
  - `NavierStokes.AxisContraction.angularLinearCoefficient`
  - `NavierStokes.AxisContraction.naturalRemainder`
  - literal branch `lin1 := O.j2 (O.product (angularLinearCoefficient O d) phi) + O.dot2 d.wStar phi + O.param2 phi d.hStar`
- `NavierStokes/AxisOperators.lean`
  - `product`
  - `regularInverse` / `j2`
  - `inverseDotProduct` / `dot2`
  - `inverseParamProduct` / `param2`

For output radial row `n>0`, all three regular-inverse branches use the pinned divisor `n(n+1)` after reading radial total degree `n-1`. `dot2` contributes the radial Euler factor carried by the second argument; `param2` differentiates the first argument once in the parameter before the inverse. Output row zero is literal zero for each branch.

## Actual-data dependency

The production factory accepts only the landed `ActualScheduleWideFirstPicardState`. From that same source-bound object it consumes:

- genuine mixed-scale `phi1` coefficient jets;
- actual SchedulePressure `AxisData.wStar`, `one`, `eta`, `uStar`, `hStar`, and scalar `h`;
- the landed actual-schedule `coefficientOperators.product` semantics for `eta*uStar`;
- the theorem-selected `Lambda` and common coefficient epsilon.

No caller-supplied coefficient table, Lambda, cutoff, derivative table, surrogate state, fitted sample, or missing-row value is admitted. The three `phi1` scale factors are never added into one binary64 value before the pinned operations.

## Regression obligations

`tests/test_axis_coefficient_wide_first_picard_lin1.py` checks:

1. identity binding to the genuine x1, actual AxisData, common epsilon, and theorem-selected Lambda;
2. exact row-zero vanishing of all three regular-inverse branches;
3. the first radial-row identity, where the `dot2` Euler factor vanishes and only `j2` plus `param2` remain;
4. a general `(n,m)` coefficient against an independently replayed radial convolution, eta-Leibniz, dot2 Euler factor, param2 eta shift, and exact `n(n+1)` divisor;
5. preservation of exactly the three genuine ordinary x1 angular scales, without inventing a pressure family;
6. rejection of surrogate states and invalid coefficient/window coordinates; and
7. machine-readable provenance with `slow1(x1)`, `lin2(x1)`, complete `naturalRemainder(x1)`, x2, fixed-point, and paper-exact claims still false.

## Remaining boundary

The genuine x1 source, pressure, complete slow2, corrected quad1, and now lin1 constituents are available on compatible actual-schedule data. The shortest remaining closure path is to materialize pinned `slow1(x1)` and `lin2(x1)`, then assemble the complete angular/axial `naturalRemainder(x1)` and only then apply the outer Picard update to produce genuine `x2`. Iteration/contraction closure, global coefficient-space membership, fixed-point `phi/u`, derived average/pressure fields, `NaturalProfileAssembly`, and paper-exact velocity remain unresolved.
