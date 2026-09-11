# Axis coefficient multiply-Y provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisWeightEstimates.lean`
  - `weight_radial_shift`
- `NavierStokes/AxisOperators.lean`
  - `rowData`
  - `multiplyYScale`
  - `multiplyYData`
  - `mulY`
  - `norm_mulY_le`
- `NavierStokes/AxisEvaluationAlgebra.lean`
  - `jet_mulY_eval`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators`
  - downstream `naturalRemainder`

## Executed identity

For compatible radial coefficient jets `J[n,m](eta)`, the pinned `multiplyYData` row map is

`J_mulY[0,m](eta) = 0`,

and, for `n >= 1`,

`J_mulY[n,m](eta) = J[n-1,m](eta)`.

This is exactly the coefficient action of multiplication by the radial variable,

`F(Y,eta) -> Y * F(Y,eta)`.

The implementation accepts only an existing `AxisCoefficientJetState`. Callers cannot choose a radial scale, epsilon, row shift, or replacement coefficient table; the zero row and predecessor-row shift are fixed by the pinned `multiplyYScale` / `multiplyYData` definitions.

## Independent checks

`tests/test_axis_coefficient_multiply_y.py` compares the production row map against direct multiplication of a finite actual-reference radial polynomial by the physical radial variable `Y`. It separately checks the first eta derivative with a centered finite-difference oracle and verifies that the actual axial reference degree-one row is shifted to degree two with unit coefficient.

## Remaining boundary

This increment materializes another genuine linear primitive required by `coefficientOperators`, but it does **not** certify the global all-index weighted norm needed for actual `AxisSpace` membership. The regular inverses, derivative/bilinear composites and any other still-missing pinned operators, the complete `coefficientOperators` record, `naturalRemainder(x0)`, the first genuine Picard iterate, fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification remain unresolved. `paper_exact_velocity_available` must remain false.
