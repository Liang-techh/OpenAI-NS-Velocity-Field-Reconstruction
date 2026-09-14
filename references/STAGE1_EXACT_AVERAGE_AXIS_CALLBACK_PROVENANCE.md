# Stage-1 exact-average axis callback precedence

## Scope

This is one downstream-only Issue #1 / Stage-1 interface increment replayed from exact repository
`main` `065b67e9a3e8e22d697b49ccc67163db60a79783` after merged PR #319 advanced the baseline.

The production and focused regression payload are unchanged from stale PR #356. The intervening
`1bc98e6... -> 065b67e...` main delta is merged Stage-2 hierarchy-owned third-eta lower-source work
and is file-disjoint from this Stage-1 profile-interface increment.

It does not implement `AxisCoefficientSpace`, coefficient operations, `naturalRemainder`, `x2`, a
fixed point, or any `phi/u/average/pressure` field. Agent 6 retains ownership of that upstream lane.

## Identity seam

`LeadingProfile` already permits an exact radial-average callback through `average_U` (and the
parallel `average_dU_deta` callback). Before this increment, `_average(...)` handled `X == 0` before
consulting that exact callback. Consequently a natural-profile adapter that supplied the backend-owned
average field could still be forced to call the underlying `U(0, eta)` at the axis.

This increment gives the supplied exact callback precedence for every admissible `X >= 0`, including
`X == 0`. The axis identity `average(U)(0, eta) = U(0, eta)` remains the fallback only when no exact
callback is supplied. Numerical Gauss-Legendre quadrature remains only the non-axis fallback when no
exact callback is available.

This is the minimal interface condition required before the separate Stage-1 exact-average handoff
can safely preserve a genuine backend-owned average field through the axis. It does not itself supply
that field.

## Verification

`tests/test_profiles_exact_average.py` contains two independent regressions:

- with `average_U` supplied, `U` is a fail-fast callable and `radial_average_U(0, eta)` must still
  return the exact callback value;
- without `average_U`, the axis path must retain the existing exact point-value fallback and call
  `U(0, eta)` exactly once.

No sampled fit, tolerance-based theorem identity, default-zero coefficient, or toy profile is used as
paper evidence.

Fresh exact-head GitHub Actions are authoritative for focused/full pytest, both Python/NumPy lanes,
wheel/outside-checkout CLI/audit/demo, the four slices, and the expected
`audit --require-paper-exact` exit 2. No CI success is inferred from the stale-base run.

## Truth boundary

Stage 1 remains `formal-structure`. Complete `AxisCoefficientSpace` closure, genuine
`naturalRemainder(x1)`, `x2`, the fixed point, materialized `phi/u/average/pressure`,
`NaturalProfileAssembly`, regular-inner-core / heat-exterior matching, and
support/moment/matching/cone closure remain open.

`full_reconstruction=false` and `paper_exact_velocity_available=false` remain mandatory.
