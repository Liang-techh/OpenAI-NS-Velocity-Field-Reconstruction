# Section 10 fixed-cutoff binding provenance

Status: **formal-structure only**. This increment does not make the reconstructed velocity paper-exact, and `paper_exact_velocity_available` must remain `false`.

## Baseline

- Integration base: `main` at `b0b8f784c66cd089436f8a1d96405ffe75b7621a`.
- Formal source: repository-pinned OpenAI `NavierStokes/SpatialLocalization.lean`.
- Pinned geometry already present in `spatial_localization.py`: `cutoffProfile(r2,z)=chi(16 r2) chi(4 z)`, plateau `r^2 < 1/32, |z| < 1/8`, closed support `r^2 <= 1/16, |z| <= 1/4`.

## What this increment adds

`section10_localized_field(local)` atomically binds the existing Section 10 executable spatial-cutoff representative to its matching analytic Cartesian gradient. Downstream callers therefore no longer need to pair the fixed cutoff with a gradient manually when using `LocalizedField`'s product rule

`curl(c A) = c curl(A) + grad(c) cross A`.

The adapter fails closed on non-`LocalField` inputs and preserves the existing zero-before-local-evaluation behavior outside the official support cylinder.

## Independent checks

The regression suite checks that the adapter selects the exact matching cutoff/gradient pair, and compares the resulting analytic product-rule velocity in the transition collar against an independently finite-differenced Cartesian curl of the localized vector potential. Exterior tests use deliberately unavailable local data and verify that points outside the fixed support return zero without evaluating that data.

These checks are independent of the analytic cutoff-gradient formula at the final curl-comparison stage; they do not set forcing equal to a residual or otherwise obtain closure tautologically.

## Non-claims

- The executable C-infinity transition collar is geometry-faithful to the pinned Mathlib `ContDiffBump` support/plateau specification, but it is not claimed to equal that noncomputable bump pointwise.
- The supplied `LocalField` is not certified here to be the completed Sections 4–9 paper field.
- No new time-localization, endpoint-limit, smooth-forcing, finite-energy, residual-convergence, or blow-up theorem is claimed.
- This adapter does not alter the repository's paper-exact gate; `paper_exact_velocity_available=false` remains required until the genuine upstream construction and Section 10 closure are complete.
