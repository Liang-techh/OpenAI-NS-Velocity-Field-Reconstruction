# Section 5 preceding-diffusion parameter-jet provenance

## Scope

The exact PositiveAxis strict-lower source already evaluates the paper/Lean quantity

`precedingDiffusion = Z_(b-D) (Z_b F_(n-1))`

from a second `(X, eta)` jet.  The next solver step toward a hierarchy-owned `partial_eta f_n` requires its eta derivative.  This increment adds only that derivative layer.

## Paper-derived structure

Differentiating `Z_(b-D)(Z_b F)` once in `eta` requires the existing second jet plus exactly three additional mixed derivatives:

- `F_XXeta`,
- `F_Xetaeta`, and
- `F_etaetaeta`.

`background_preceding_diffusion_parameter_jet.py` exposes these as `ProfileThirdMixedJet` and returns both the existing preceding-diffusion value and its analytic eta derivative.  Production performs the quotient/product differentiation explicitly; it does not finite-difference eta and does not fit a sampled coefficient.

The value is delegated back to the already-landed `preceding_diffusion_from_second_jet` implementation so the new derivative layer cannot silently fork the exact value path.

## Verification

`tests/test_background_preceding_diffusion_parameter_jet.py` uses a bivariate polynomial of degree two in `X` and three in `eta`, so all required third-mixed derivatives are nontrivial.  For several positive orders, both angular/axial paper powers, and points including `X=0`, the analytic eta derivative is compared against a centered finite difference of the older value-only production path with independently recomputed exact polynomial jets.

A separate regression perturbs `F_XXeta`, `F_Xetaeta`, and `F_etaetaeta` one at a time and verifies that every newly required derivative affects the result.  Missing stronger jets and invalid domain/order data fail closed.

## Truth boundary

This remains Stage-2 `formal-structure` infrastructure and `paper_exact_velocity_available=false` remains mandatory.

The real hierarchy does not yet own `ProfileThirdMixedJet` for `phi_(n-1)` through the Lemma 5.2 compact repair, and Issue #1 still owns the real leading profile.  Therefore this increment does **not** claim hierarchy-owned `partial_eta actualLowerSource`, hierarchy-owned `partial_eta f_n`, the genuine `k=1` Eq. (5.7) Picard application, Picard convergence, a recursively materialized all-order coefficient family, the paper's recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity.
