# Axis coefficient wide first-Picard average provenance

## Scope

This increment materializes exactly one new operator seam needed to evaluate the pinned natural-axis fixed-point remainder away from the reference center: the radial average of the genuine theorem-scale first Picard axial field,

`b u1 = coefficientOperators.average(u1)`.

The source theorem structure is pinned to `openai/NavierStokesAndEuler` commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, specifically the `AxisOperators.averageData` row action and the `AxisContraction.naturalRemainder` assignment `bu := O.average u`.

## Exact executable identity

For every radial degree `n >= 0`, parameter-jet order `m`, and admissible `eta`, the pinned average is row-local:

`J_avg[n,m](eta) = J[n,m](eta) / (n + 1)`.

The already-landed first Picard axial coefficient is stored losslessly as four pieces:

`u1 = u0 + B/(2 Lambda) + S/(2 Lambda^2) + P/(2 Lambda)`,

where the pressure piece is carried in split signed-log form. The lifted average therefore preserves the decomposition exactly:

`average(u1) = average(u0) + average(B)/(2 Lambda) + average(S)/(2 Lambda^2) + average(P)/(2 Lambda)`.

At one row this is implemented by dividing each ordinary Decimal numerator by the exact integer `n+1`. For a nonzero signed-log pressure jet, the sign and common amplitude log are unchanged and `log(n+1)` is subtracted only from the separate factor log. Exact zero remains exact zero.

No binary64 projection of the theorem-scale correction is required and no caller may supply a replacement radial scale, Lambda, pressure term, coefficient table, derivative table, or surrogate state.

## Independent regression boundary

The tests:

- require the state to be bound to the landed genuine `ActualScheduleWideFirstPicardState`;
- check the exact `1/(n+1)` Decimal action on the reference and both inverse-Lambda numerators;
- cross-check the O(1) reference component against the already-landed independent `axis_coefficient_average` implementation;
- verify the signed-log pressure action by the literal `-log(n+1)` identity while preserving sign/common scale;
- require nonzero out-of-binary64 pressure to fail closed rather than collapse to zero;
- verify row zero is the exact identity factor and preserves exact zero pressure;
- reject negative radial degree and non-first-Picard input types.

## Truth boundary

This is **not** `naturalRemainder(x1)` and is not `x2` or a fixed-point certificate. It materializes only the first operator assignment `bu = average(u)` for the mixed-scale `u1` input. The mixed-scale product, inverse-bilinear, slow, source/pressure, resolvent, and final recombination operations required for the complete nonlinear remainder at `x1` are still unresolved.

The global all-index weighted `AxisSpace` membership/norm certificate is also not supplied here. The underlying actual-schedule phase still inherits the landed numerical quadrature boundary. Therefore Stage 1 remains `formal-structure`, `paper_exact_velocity_available=false`, and `full_reconstruction=false`.
