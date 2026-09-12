# Axis coefficient naturalRemainder provenance

## Source pin

This increment is derived from `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, module
`NavierStokes/AxisContraction.lean`, specifically the definitions
`angularLinearCoefficient`, `angularQuadraticCoefficient`, `averageCoefficient`,
`angularSlowCoefficient`, `axialLinearCoefficient`,
`axialQuadraticCoefficient`, and `naturalRemainder`.

## Executable construction

`src/openai_ns_reconstruction/axis_coefficient_natural_remainder.py` composes
only already-landed actual-scale objects:

- `ActualScheduleAxisCoefficientData` for the fixed `A,D,h` and coefficient fields;
- the complete pinned `AxisCoefficientOperators` record;
- the coefficientwise exact `AxisCoefficientNaturalResolvent`;
- lazy `AxisCoefficientJetState` inputs for the theorem variables `a`, `phi`, and `u`.

The code follows the Lean `let` structure term-by-term: `bu`, `lin1`, `quad1`,
`slow1`, `lin2`, `slow2`, `source`, and `pressure`, followed by the angular
resolvent and the axial `inverseL` multiplication.  No operator table,
pressure datum, sigma, epsilon, or fixed coefficient field can be replaced by
a caller.

## Validation

The regression suite checks consequences of the pinned expression that are
independent of any paper-exact amplitude choice:

1. radial row zero vanishes because every integrated remainder term passes
   through the pinned regular inverse;
2. the angular component is independent of `a`, while the axial pressure
   contribution is quadratic in `a` because `source=(a*a)*(phi*phi)`;
3. both components are affine in `t`, since `t` occurs only in `-t*slow1` and
   `-t*slow2`;
4. coefficient states with a mismatched epsilon are rejected.

The zero/scaled amplitude states used in algebraic regression are explicitly
marked test-only.  They are not manuscript data and are never promoted to a
paper-exact construction.

## Remaining boundary

This increment materializes the pinned **composition**, but not yet the final
manuscript instantiation.  In particular it does not yet construct the
coefficient state of `realAmplitude h j sigma Lambda C`, does not set
`t=1/Lambda` from the theorem-selected wide scale chain, and therefore does not
yet evaluate `naturalRemainder(x0)` or the first genuine Picard iterate.
Global all-index weighted `AxisSpace` membership/norm certification also remains
open.  Consequently `paper_exact_velocity_available=false` remains mandatory.
