# Axis coefficient inverse parameter-product provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisOperators.lean`
  - `differentialFamily`
  - `inverseParamProductData`
  - `inverseParamProduct`
  - `norm_inverseParamProduct_le`
  - `jet_inverseParamProduct`
- `NavierStokes/AxisEvaluationAlgebra.lean`
  - `inverseParamProduct_equation`
- `NavierStokes/NaturalAxisBridge.lean`
  - evaluated identity identifying the operator with `partialEta(profile A) * profile B`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators` entries `param1` and `param2`
  - downstream `naturalRemainder`

## Executed identity

For compatible coefficient jets `J_A[n,m](eta)`, `J_B[n,m](eta)` and positive integer `r`, the pinned specialization `differentialFamily r 1 (fun _ => 1)` gives

`J_out[0,m](eta) = 0`,

and, for `n >= 1`,

`J_out[n,m](eta) = (1 / radialDivisor(r,n-1)) *`
`  sum_{i+j=n-1} sum_{k+l=m} binom(m,k)`
`  J_A[i,k+1](eta) J_B[j,l](eta)`.

Thus the evaluated operator is the zero-datum regular radial inverse `J_r` applied to `(partial_eta A) * B`. The two instances used by the pinned `AxisContraction.coefficientOperators` record are `r=1` (`param1`) and `r=2` (`param2`).

The implementation accepts only two existing `AxisCoefficientJetState` operands and the theorem-domain positive integer `r`. It requires the exact same pinned epsilon on both inputs. Callers cannot inject a parameter-derivative table, radial divisor table, epsilon, integration constant, or replacement coefficient table.

## Independent checks

`tests/test_axis_coefficient_inverse_param_product.py` exercises the actual SchedulePressure-derived reference states. Its main cross-check constructs an exact parameter-derivative view from the landed actual jets, passes that through the separately landed `productFamily` implementation and then the separately landed `regularInverse`, and compares the result against the new direct `differentialFamily` specialization for `r=1,2`, several radial rows, and several eta-jet orders. A second actual-reference regression uses the fact that the axial reference has only radial degree one to check closed low-order formulas. The tests also differentiate the zeroth output coefficient in eta by centered finite difference and independently check the physical radial equation using finite-differenced actual axial coefficients rather than the production `m=1` source jets. Invalid weight-scale mixtures, nonpositive/noninteger `r`, invalid indices, and eta values outside the pinned window fail closed.

## Remaining boundary

This increment materializes the `inverseParamProduct` primitive required by `coefficientOperators`, including the `param1` and `param2` actions, but it does **not** certify the global all-index weighted norm needed for actual `AxisSpace` membership. The inverse dot-product and inverse mixed bilinear composites, the complete `coefficientOperators` record, `naturalRemainder(x0)`, the first genuine Picard iterate, fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification remain unresolved. `paper_exact_velocity_available` must remain false.
