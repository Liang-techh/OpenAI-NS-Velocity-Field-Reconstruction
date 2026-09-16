# Axis coefficient wide first-Picard source provenance

## Scope

This increment materializes only the pinned source constituent
`source(x1) = a^2 * phi1^2` for the genuine actual-SchedulePressure first Picard
state. It consumes the already-landed typed `ActualScheduleWideFirstPicardState`,
the landed mixed-scale `phi1 * phi1` convolution, and the theorem-selected
signed-log amplitude carried by the same x1 chain.

## Pinned source

Formal source repository: `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant definitions are `NavierStokes.AxisContraction.naturalRemainder` and the
coefficient product in `NavierStokes.AxisOperators`. The natural remainder
source is `a^2 * phi^2`. At x1 the landed angular field is already represented
as `phi1 = r + b/Lambda + c/Lambda^2`; its exact coefficient product therefore
has powers Lambda^0 through Lambda^-4.

## Executable realization

`axis_coefficient_wide_first_picard_source.py` keeps the common `a^2` scale out
of binary64. Since the amplitude has radial degree zero, each radial row of the
landed `phi1^2` state is preserved. Eta derivatives use the exact finite
Leibniz sum. The normalized amplitude derivatives `d^k(a^2)/a^2` are generated
by the complete Bell recurrence from the same actual normalized-gradient jets
and the same wide Decimal Lambda already bound into x1.

The output contains five Decimal factors, one for each inverse-Lambda power
0..4, plus the exact common amplitude log. `source_terms_log()` restores each
term as `a^2 * s_p / Lambda^p` in split signed-log coordinates. No caller may
supply a coefficient table, pressure row, Lambda, C, amplitude, derivative
sample, or replacement x1.

## Regression boundary

The regression independently checks the zeroth-eta identity against the landed
`phi1^2` state and the first-eta product rule
`(a^2 phi1^2)'/a^2 = (phi1^2)' + 2 Lambda realGradient * phi1^2` at explicit
96-digit Decimal precision. It also checks that all nonzero source powers keep
`2*log(a)` as the common signed-log scale and the correct Lambda denominator.

This does **not** establish the x1 pressure chain, angular/axial recombination,
complete `naturalRemainder(x1)`, `x2`, Picard convergence, global
`AxisCoefficientSpace` membership, a fixed point, `NaturalProfileAssembly`, or
paper-exact velocity. The underlying phase evaluation remains the previously
documented numerical-quadrature boundary.
