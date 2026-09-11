# Axis coefficient regular-inverse provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisWeightEstimates.lean`
  - `radialDivisor`
  - `regularInverseJet`
  - `regularInverseJet_bound`
- `NavierStokes/AxisOperators.lean`
  - `inverseScale`
  - `inverseData`
  - `regularInverse`
  - `norm_regularInverse_le`
  - `jet_regularInverse_zero`
  - `jet_regularInverse_succ`
- `NavierStokes/AxisEvaluationAlgebra.lean`
  - `regularInverse_equation`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators` entries `j1` and `j2`
  - downstream `naturalRemainder`

## Executed identity

For compatible radial coefficient jets `J[n,m](eta)` and a positive integer `r`, the pinned regular zero-datum inverse has

`J_inv[0,m](eta) = 0`,

and, for `n >= 1`,

`J_inv[n,m](eta) = J[n-1,m](eta) / radialDivisor(r,n-1)`,

with

`radialDivisor(r,k) = (k+1)(k+r)`.

At the evaluated-profile level this is the regular zero-datum right inverse of

`Y g'' + r g' = f`.

The implementation accepts only an existing `AxisCoefficientJetState` and the theorem-domain positive integer `r`. Callers cannot inject a divisor table, integration constant, epsilon, radial scale, or replacement coefficient table. The two instances needed by the pinned `AxisContraction.coefficientOperators` record are `r=1` (`j1`) and `r=2` (`j2`).

## Independent checks

`tests/test_axis_coefficient_regular_inverse.py` compares the production coefficient map for `r=1` and `r=2` against an independent Green-kernel quadrature of the zero-datum radial ODE on finite actual-reference radial polynomials. It separately checks the physical equation `Y g'' + r g' = f`, an eta derivative by centered finite difference, and the exact degree-one to degree-two scaling of the actual axial reference state. Invalid nonpositive/noninteger `r`, invalid indices, and out-of-window eta values fail closed.

## Remaining boundary

This increment materializes the pinned regular-inverse primitive required by `coefficientOperators`, including the `j1` and `j2` row actions, but it does **not** certify the global all-index weighted norm needed for actual `AxisSpace` membership. The inverse parameter-product / inverse dot-product / inverse mixed bilinear composites and any other still-missing pinned operators, the complete `coefficientOperators` record, `naturalRemainder(x0)`, the first genuine Picard iterate, fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification remain unresolved. `paper_exact_velocity_available` must remain false.
