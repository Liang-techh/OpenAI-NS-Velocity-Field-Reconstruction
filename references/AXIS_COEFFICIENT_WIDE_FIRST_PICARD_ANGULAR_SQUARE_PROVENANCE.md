# Axis coefficient wide first-Picard angular-square provenance

## Scope

This increment materializes exactly one nonlinear product seam needed to evaluate the pinned natural-axis remainder at the genuine theorem-scale first Picard iterate:

`phi1 * phi1`.

The source algebra is pinned to `openai/NavierStokesAndEuler` commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, specifically `AxisWeightEstimates.jetProduct` / `AxisOperators.productFamily` and the `AxisContraction.naturalRemainder` source expression `a * a * phi * phi`.

## Exact executable identity

The landed first Picard angular jet is kept as

`phi1 = r + b/Lambda + c/Lambda^2`.

For radial degree `n`, eta-derivative order `m`, and admissible `eta`, the pinned product law is

`J_{AB}[n,m] = sum_{i+j=n} sum_{k+l=m} binom(m,k) J_A[i,k] J_B[j,l]`.

Applying that finite identity to both copies of `phi1` produces powers from `Lambda^0` through `Lambda^-4`. The implementation collects each power in an independent 96-digit Decimal numerator and never adds the tiny inverse-Lambda corrections into the O(1) reference component. Thus a mathematically nonzero theorem-scale correction cannot disappear merely because `Lambda` is far outside binary64 scale.

No caller may provide a replacement Lambda, angular coefficient table, surrogate first-Picard state, fitted derivative data, or product cutoff. The state accepts only the already-landed genuine `ActualScheduleWideFirstPicardState`, which itself is bound to the actual SchedulePressure theorem-selection chain.

## Independent regression boundary

The tests:

- require the product state to be bound to the genuine theorem-scale `x1`;
- verify the row-zero product against the closed five-term expansion `r^2`, `2rb`, `b^2+2rc`, `2bc`, `c^2`;
- cross-check the complete O(1) radial/eta convolution against the independently landed binary64 `axis_coefficient_product(reference.phi, reference.phi)` implementation;
- independently reconstruct every inverse-Lambda numerator from the literal radial convolution and eta-Leibniz sum;
- verify the correction terms remain separate after division by `Lambda`, `Lambda^2`, `Lambda^3`, and `Lambda^4`;
- reject invalid indices and non-first-Picard inputs.

## Truth boundary

This is **not** the complete `naturalRemainder(x1)`. It materializes only one nonlinear product needed by the next pressure/source evaluation and establishes the scale bookkeeping required for later mixed-scale bilinear operators. Products involving the axial first-Picard field still require signed-log pressure-aware arithmetic, and the remaining linear/parameter/mixed/inverse/resolvent operations have not yet been lifted to `x1`.

No Picard `x2`, fixed-point convergence certificate, final `phi/u`, derived average/pressure, `NaturalProfileAssembly`, support/moment/matching/cone closure, or global all-index weighted `AxisSpace` certificate is claimed. The underlying actual-schedule phase still inherits the landed numerical quadrature boundary. Therefore Stage 1 remains `formal-structure`, `paper_exact_velocity_available=false`, and `full_reconstruction=false`.
