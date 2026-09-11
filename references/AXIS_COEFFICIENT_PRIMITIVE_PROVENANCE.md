# Axis coefficient radial-primitive provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisOperators.lean`
  - `rowData`
  - `primitiveScale`
  - `primitiveData`
  - `primitive`
  - `norm_primitive_le`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators`
  - downstream `naturalRemainder`

## Executed identity

For compatible radial coefficient jets `J[n,m](eta)`, the pinned `primitiveData` row map is

`J_prim[0,m](eta) = 0`,

and, for `n >= 1`,

`J_prim[n,m](eta) = J[n-1,m](eta) / n`.

This is the coefficient identity for the zero-at-axis radial primitive

`P[F](Y,eta) = integral_0^Y F(s,eta) ds`.

The implementation accepts only an existing `AxisCoefficientJetState`. Callers cannot choose a replacement integration constant, epsilon, row scale, or coefficient table; the zero constant is fixed by the pinned `primitiveScale 0 = 0` row.

## Independent checks

`tests/test_axis_coefficient_primitive.py` compares the production row map against SciPy adaptive quadrature of a finite radial polynomial built from the actual SchedulePressure-derived angular reference coefficients. It separately checks a parameter derivative with a centered finite-difference oracle and verifies that the actual axial reference degree-one row is moved to degree two with the exact factor `1/2`.

## Remaining boundary

This increment materializes another genuine linear primitive required by `coefficientOperators`, but it does **not** certify the global all-index weighted norm needed for actual `AxisSpace` membership. `parameterPrimitive`, the regular inverses / derivative composites and any other still-missing pinned operators, the complete `coefficientOperators` record, `naturalRemainder(x0)`, the first genuine Picard iterate, fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification remain unresolved.
