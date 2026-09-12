# Axis coefficient inverse mixed provenance

Status: **formal-structure only**. This artifact does not make the Stage 1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisOperators.lean`
  - `differentialFamily`
  - `inverseMixedData`
  - `inverseMixed`
  - `norm_inverseMixed_le`
  - `jet_inverseMixed`
- `NavierStokes/AxisEvaluationAlgebra.lean`
  - evaluated inverse-mixed equation
- `NavierStokes/NaturalAxisBridge.lean`
  - evaluated identity identifying the operator with `(partial_eta profile A) * (Y * partialY(profile B))`
- `NavierStokes/AxisContraction.lean`
  - `coefficientOperators` entries `mixed1` and `mixed2`
  - downstream `naturalRemainder`

## Executed identity

For compatible coefficient jets `J_A[n,m](eta)`, `J_B[n,m](eta)` and positive integer `r`, the pinned specialization `differentialFamily r 1 (fun j => j)` gives

`J_out[0,m](eta) = 0`,

and, for `n >= 1`,

`J_out[n,m](eta) = (1 / radialDivisor(r,n-1)) *`
`  sum_{i+j=n-1} sum_{k+l=m} j * binom(m,k)`
`  J_A[i,k+1](eta) J_B[j,l](eta)`.

The `k+1` shift is exactly the parameter derivative on the first operand, while the factor `j` is exactly the coefficient action of the Euler radial derivative `Y * partial_Y` on the second operand. Thus the evaluated operator is the zero-datum regular radial inverse `J_r` applied to `(partial_eta A) * (Y * partial_Y B)`. The two instances used by the pinned `AxisContraction.coefficientOperators` record are `r=1` (`mixed1`) and `r=2` (`mixed2`). The pinned Lean norm estimate is `||inverseMixed|| <= 5120 / epsilon`; this Python increment materializes the exact jet action but does not claim the all-index weighted norm certificate required for genuine `AxisSpace` membership.

The implementation accepts only two existing `AxisCoefficientJetState` operands and the theorem-domain positive integer `r`. It requires the exact same pinned epsilon on both inputs. Callers cannot inject a parameter/radial derivative table, radial divisor table, epsilon, integration constant, or replacement coefficient table.

## Independent checks

`tests/test_axis_coefficient_inverse_mixed.py` exercises the actual SchedulePressure-derived reference states. Its main cross-check constructs two test-side exact views: `partial_eta A` by shifting the landed parameter-jet index, and `Y * partial_Y B` by multiplying radial row `n` by `n`. Those views are passed through the separately landed `productFamily` implementation and then the separately landed `regularInverse`, and the result is compared against the new direct `differentialFamily(p=1,d(j)=j)` specialization for `r=1,2`, several radial rows, and several eta-jet orders. A second actual-reference regression uses the fact that the axial reference has only radial degree one to check closed low-order formulas. The tests also differentiate the zeroth output coefficient in eta by centered finite difference and independently reconstruct `(partial_eta u)(Y) * (Y * partial_Y phi)(Y)` to check the physical radial equation `Y g'' + r g' = source`. Invalid weight-scale mixtures, nonpositive/noninteger `r`, invalid indices, and eta values outside the pinned window fail closed.

## Remaining boundary

This increment materializes the last previously missing inverse-bilinear primitive required by the pinned `coefficientOperators` record, including the `mixed1` and `mixed2` actions. It does **not** yet assemble and validate that complete record as one executable object and does not certify the global all-index weighted norm needed for actual `AxisSpace` membership. The next direct Stage 1 blocker is therefore assembly of the pinned `coefficientOperators` interface from the landed product/average/primitive/parameterPrimitive/mulY/j1/j2/param1/param2/dot1/dot2/mixed1/mixed2 actions, followed by the first genuine evaluation of `naturalRemainder(x0)`. The real Picard iterate, fixed-point `phi/u`, average/pressure reconstruction, `NaturalProfileAssembly`, and final support/moment/matching/cone verification remain unresolved. `paper_exact_velocity_available` must remain false.
