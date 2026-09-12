# Axis coefficient inverse dot-product provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisOperators.lean`
  - `differentialFamily`
  - `inverseDotProductData`
  - `inverseDotProduct`
  - `norm_inverseDotProduct_le`
  - `jet_inverseDotProduct`
- `NavierStokes/AxisEvaluationAlgebra.lean`
  - `inverseDotProduct_equation`
- `NavierStokes/NaturalAxisBridge.lean`
  - evaluated identity identifying the operator with `profile A * (Y * partialY(profile B))`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators` entries `dot1` and `dot2`
  - downstream `naturalRemainder`

## Executed identity

For compatible coefficient jets `J_A[n,m](eta)`, `J_B[n,m](eta)` and positive integer `r`, the pinned specialization `differentialFamily r 0 (fun j => j)` gives

`J_out[0,m](eta) = 0`,

and, for `n >= 1`,

`J_out[n,m](eta) = (1 / radialDivisor(r,n-1)) *`
`  sum_{i+j=n-1} sum_{k+l=m} j * binom(m,k)`
`  J_A[i,k](eta) J_B[j,l](eta)`.

The factor `j` is exactly the coefficient action of the Euler radial derivative `Y * partial_Y` on the second operand. Thus the evaluated operator is the zero-datum regular radial inverse `J_r` applied to `A * (Y * partial_Y B)`. The two instances used by the pinned `AxisContraction.coefficientOperators` record are `r=1` (`dot1`) and `r=2` (`dot2`).

The implementation accepts only two existing `AxisCoefficientJetState` operands and the theorem-domain positive integer `r`. It requires the exact same pinned epsilon on both inputs. Callers cannot inject a radial-derivative table, radial divisor table, epsilon, integration constant, or replacement coefficient table.

## Independent checks

`tests/test_axis_coefficient_inverse_dot_product.py` exercises the actual SchedulePressure-derived reference states. Its main cross-check constructs an exact test-side coefficient view of `Y * partial_Y B` by multiplying radial row `n` by `n`, passes that view through the separately landed `productFamily` implementation and then the separately landed `regularInverse`, and compares the result against the new direct `differentialFamily` specialization for `r=1,2`, several radial rows, and several eta-jet orders. A second actual-reference regression uses the fact that the axial reference has only radial degree one to check closed low-order formulas. The tests also differentiate the zeroth output coefficient in eta by centered finite difference and independently reconstruct `phi(Y) * (Y * partial_Y u(Y))` to check the physical radial equation `Y g'' + r g' = source`. Invalid weight-scale mixtures, nonpositive/noninteger `r`, invalid indices, and eta values outside the pinned window fail closed.

## Remaining boundary

This increment materializes the `inverseDotProduct` primitive required by `coefficientOperators`, including the `dot1` and `dot2` actions, but it does **not** certify the global all-index weighted norm needed for actual `AxisSpace` membership. The inverse mixed bilinear composite, the complete `coefficientOperators` record, `naturalRemainder(x0)`, the first genuine Picard iterate, fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification remain unresolved. `paper_exact_velocity_available` must remain false.
