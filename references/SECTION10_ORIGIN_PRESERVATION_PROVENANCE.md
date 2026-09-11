# Section 10 late-origin / blow-up-preservation provenance

Status: **formal-structure**.

This increment records only the localization implication that preserves an
already-certified incoming origin blow-up.  It does not construct, sample-fit,
or certify the unresolved paper velocity.

## Paper location

- Section 10 final localization, especially Eq. (10.4).
- The blow-up path recorded in the repository source ledger at Eqs. (10.20)-(10.21).

## Official Lean cross-check

Pinned source: `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, module
`NavierStokes/SpatialLocalization.lean`.

The relevant theorem chain is:

- `zero_mem_plateau` and `spatialCutoff_eq_one`: the origin lies in the fixed
  spatial plateau, so the cutoff is locally exactly one there;
- `periodicVelocity_eventuallyEq_cut` / `periodicVelocity_eq`: in the central
  no-overlap cube, periodization agrees locally with the cut field;
- `localizedVelocity_eq`: for `t >= 3/4` and spatial points in the plateau, the
  final localized velocity equals the incoming `SpatialCurl.spatialCurl A`;
- `localizedVelocity_origin`: the preceding identity specialized to the origin;
- `localizedVelocity_origin_blowup`: a proved left-limit
  `||spatialCurl A(t,0)|| -> infinity` as `t -> 1-` is preserved by final
  localization;
- `localizedVelocity_speed_unbounded`: that origin limit implies the formal
  `SpeedUnboundedAtOne` conclusion.

The time factor is the already-landed `TimeLocalization.timeSwitch`: it is
exactly one for the late branch `t >= 3/4`.

## Executable artifact

`src/openai_ns_reconstruction/section10_origin_preservation.py` adds a
fail-closed `Section10LateOriginCertificate` for `3/4 <= t < 1`.  It checks,
using the current executable Section 10 modules, that at the origin:

- the point is in the official plateau geometry;
- the spatial cutoff value is exactly `1`;
- the analytic spatial cutoff gradient is exactly zero;
- the time switch is exactly `1`.

It also evaluates the cutoff product-rule value

`c curl(A) + grad(c) x A`

for caller-supplied finite point data and therefore reduces exactly to
`curl(A)` at the certified late origin.  A separate helper records that the
late time activation also leaves an already-known origin value unchanged.

`tests/test_section10_origin_preservation.py` independently reads the landed
spatial/time localization functions, checks the exact factors at the late
endpoint branch (including the last representable float below `t=1`), checks
that arbitrary finite algebra data are unchanged by the product-rule/time
factors, and verifies fail-closed behavior outside `3/4 <= t < 1` or for
malformed/nonfinite vectors.

## Why the executable cutoff representative is sufficient for this increment

The repository does **not** claim pointwise equality with Mathlib's
noncomputable `ContDiffBump` in the transition collars.  That unresolved collar
issue does not affect this particular identity: the origin is strictly inside
the common plateau where both specifications are exactly one, its cutoff
gradient is zero there, and the late time branch is likewise identically one.
Thus this increment tests only shared exact plateau geometry and never uses a
collar value as paper data.

## Boundary / blocker

This certificate does **not** establish the premise

`||SpatialCurl.spatialCurl A (t,0)|| -> infinity` as `t -> 1-`.

That premise must come from the completed upstream Sections 4-9 construction.
In particular, supplying arbitrary point vectors to the product-rule helper is
not evidence of the paper field, a blow-up sequence, or a divergent limit.
The module intentionally has no sampled-slope or threshold heuristic for
promoting finite values into a blow-up claim.

It also does not materialize the missing local field, prove the uniform kinetic
energy bound, produce the actual residual endpoint derivative family, obtain
the smooth global force extension through `t=1`, or independently close the
forced Navier-Stokes equation.  Those Issue #4 blockers remain unchanged.

Accordingly Stage 7 remains **formal-structure** and
`paper_exact_velocity_available` must remain **false**.
