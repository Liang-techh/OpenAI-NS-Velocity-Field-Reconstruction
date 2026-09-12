# Section 5 preceding-diffusion second parameter-jet provenance

## Scope

The PositiveAxis strict-lower source contains

`precedingDiffusion = Z_(b-D) (Z_b F_(n-1))`.

The landed first-parameter bridge already returns the exact value together with its analytic `partial_eta` derivative.  The next prerequisite for a hierarchy-owned `partial_eta^2 actualLowerSource` is the second eta derivative of this same term.  This increment adds only that stronger derivative layer.

## Paper-derived structure

Differentiating `Z_(b-D)(Z_b F)` twice in `eta` requires the existing `ProfileThirdMixedJet` plus exactly three additional eta-bearing total-order-four entries:

- `F_XXetaeta`,
- `F_Xetaetaeta`, and
- `F_etaetaetaeta`.

`background_preceding_diffusion_second_parameter_jet.py` exposes these as `ProfileFourthMixedJet`.  Its value and first eta derivative are delegated to the previously landed `preceding_diffusion_parameter_jet` after exact projection to the third-mixed jet; only the new second derivative is evaluated by the explicit analytic product/quotient rules.  Production performs no finite differencing and fits no sampled coefficient.

## Verification

`tests/test_background_preceding_diffusion_second_parameter_jet.py` uses an independent bivariate polynomial of degree two in `X` and four in `eta`, so every newly required fourth-mixed derivative is nonzero.  Across several positive recursive orders, both angular/axial paper powers, and points including `X=0`, the new bridge must reproduce the already-landed value and first derivative exactly.  Its analytic second derivative is then compared against a centered finite difference of the older analytic first-derivative path, with exact polynomial jets recomputed at the displaced eta points; finite differencing is test-only.

A separate regression perturbs `F_XXetaeta`, `F_Xetaetaeta`, and `F_etaetaetaeta` one at a time and verifies that each new entry affects the second derivative.  Weaker jets and non-finite stronger data fail closed.

## Truth boundary

This remains Stage-2 `formal-structure` infrastructure and `paper_exact_velocity_available=false` remains mandatory.

`ProfileFourthMixedJet` is an explicit stronger input contract, not a claim that the repaired hierarchy already owns this data.  In particular, caller/test polynomial jets are not paper-exact profile data.  The complete hierarchy-owned `partial_eta^2 actualLowerSource` and `partial_eta^2 f_n` still require corresponding stronger derivatives for the other strict-lower terms, including the Eq. (5.6) omega quotient path.  Issue #1 still owns the real leading profile.  Therefore this increment does **not** claim hierarchy-owned `partial_eta W_n^(1)`, further Picard iteration or convergence, recursively materialized all-order coefficients, the paper's recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity.
