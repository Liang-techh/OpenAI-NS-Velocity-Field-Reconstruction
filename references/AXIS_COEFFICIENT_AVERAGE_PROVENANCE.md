# Axis coefficient radial-average provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisOperators.lean`
  - `rowData`
  - `averageData`
  - `average`
  - `norm_average_le`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators`
  - downstream `naturalRemainder`
- `NavierStokes/NaturalAxisBridge.lean`
  - profile bridge showing `AxisOperators.average` applied to the axial coefficient state in the angular remainder.

## Executed identity

For the compatible radial coefficient jets `J[n,m](eta)`, the pinned `averageData` is

`J_avg[n,m](eta) = J[n,m](eta) / (n + 1)`.

This is the coefficient identity for the regular radial average

`A[F](Y,eta) = integral_0^1 F(t Y,eta) dt`.

The implementation accepts only an existing `AxisCoefficientJetState`; callers cannot substitute an independent epsilon, row scale, coefficient table, or fitted average.

## Independent checks

`tests/test_axis_coefficient_average.py` checks the production row map against a SciPy adaptive quadrature of a finite radial polynomial built from the actual SchedulePressure-derived reference coefficients. It also checks a parameter derivative using an independent centered finite difference and verifies the exact degree-one factor `1/2` for the actual axial reference state.

## Remaining boundary

This increment materializes one more genuine linear primitive used by `coefficientOperators`, but it does **not** certify the global all-index weighted norm needed for actual `AxisSpace` membership. The remaining pinned linear primitives/composites, the complete `coefficientOperators` record, `naturalRemainder(x0)`, the first genuine Picard iterate, the fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification are still unresolved.
