# Axis coefficient wide first-Picard axial-square provenance

## Scope

This increment materializes exactly one signed-log-aware nonlinear product seam needed to evaluate the pinned natural-axis remainder at the genuine theorem-scale first Picard iterate:

`u1 * u1`.

The source algebra is pinned to `openai/NavierStokesAndEuler` commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, specifically `AxisWeightEstimates.jetProduct` / `AxisOperators.productFamily` and the `AxisContraction.naturalRemainder.slow2` occurrence `axialQuadraticCoefficient * (u * u)`.

## Exact executable scale decomposition

The landed first Picard axial coefficient jet is already represented without destructive summation as

`u1 = r + b/Lambda + c/Lambda^2 + pressure/(2 Lambda)`.

The landed wide axial pressure has the exact coefficientwise form

`d_eta^m pressure_n(eta) = a(eta)^2 q[n,m](eta)`,

where the common theorem-selected amplitude is

`a(eta) = exp(Lambda * realPhase(eta)) / C`

and `q` is the Decimal normalized factor produced by the already-landed actual SchedulePressure pressure chain. Therefore each first-Picard pressure jet can be kept as

`a^2 p/Lambda`, with `p=q/2`,

without narrowing `a`, `1/Lambda`, or the pressure coefficient into binary64.

Applying the pinned radial convolution and eta-Leibniz identity

`J_{AB}[n,m] = sum_{i+j=n} sum_{k+l=m} binom(m,k) J_A[i,k] J_B[j,l]`

to both copies of `u1` yields the finite scale decomposition

`u1^2 = o0 + o1/Lambda + o2/Lambda^2 + o3/Lambda^3 + o4/Lambda^4`

`       + a^2 (q1/Lambda + q2/Lambda^2 + q3/Lambda^3)`

`       + a^4 s2/Lambda^2`.

The five ordinary numerators, three pressure-linear numerators, and one pressure-square numerator are accumulated independently in 96-digit Decimal arithmetic. The common `a^2` and `a^4` magnitudes are not evaluated directly. When an actual pressure term is requested, it is exposed as `SignedLogCoefficientJet` with the enormous amplitude logarithm retained in `log_scale` and only the normalized numerator plus explicit inverse-Lambda power placed in `log_factor`.

No caller can replace `Lambda`, `C`, the amplitude, pressure datum, pressure coefficients, first-Picard coefficient table, derivative table, or radial/product cutoff. The constructor accepts only the already-landed genuine `ActualScheduleWideFirstPicardState` and follows its actual SchedulePressure theorem-selection chain.

## Independent regression boundary

The tests:

- require the product state to be bound to the genuine theorem-scale `x1`;
- cross-check the complete O(1) radial/eta convolution against the independently landed binary64 `axis_coefficient_product(reference.u, reference.u)` implementation;
- independently reconstruct all five ordinary, all three pressure-linear, and the pressure-square numerators from the literal radial convolution and eta-Leibniz sum;
- verify ordinary inverse-Lambda corrections remain separate instead of being added into the O(1) reference value;
- locate an actual nonzero wide-pressure row and verify the pressure-linear terms carry the exact common `a^2` log scale, while the pressure-square term carries the exact common `a^4` log scale;
- verify the explicit `Lambda^-1`, `Lambda^-2`, and `Lambda^-3` corrections appear only in the moderate signed-log factor component;
- verify current nonzero theorem-scale pressure-product terms fail closed on binary64 underflow rather than silently becoming zero;
- reject invalid indices, out-of-window eta values, and surrogate/non-first-Picard inputs.

## Truth boundary

This is **not** the complete `naturalRemainder(x1)`. It closes the hardest first pressure-aware quadratic product used by `slow2`, but other mixed-scale products and bilinear/parameter/mixed/inverse/resolvent operations still have to be lifted before the entire remainder can be evaluated at `x1`.

No Picard `x2`, fixed-point convergence certificate, final `phi/u`, derived average/pressure, `NaturalProfileAssembly`, support/moment/matching/cone closure, or global all-index weighted `AxisSpace` certificate is claimed. The underlying actual-schedule phase still inherits the landed numerical quadrature boundary. Therefore Stage 1 remains `formal-structure`, `paper_exact_velocity_available=false`, and `full_reconstruction=false`.
