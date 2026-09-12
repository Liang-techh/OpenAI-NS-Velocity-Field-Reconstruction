# Section 5 Eq. (5.6) omega second parameter-jet provenance

## Scope

The regular Section 5 pressure/source chain uses the quotient `Omega_k/X` from Eq. (5.6).  The landed `background_omega_parameter_jet.py` already returns its value and analytic first eta derivative.  A hierarchy-owned second forcing jet will eventually require `partial_eta^2(Omega_k/X)`.  This increment adds only that stronger analytic derivative layer.

## Paper-derived structure

Writing `V_j = X beta_j`, the landed Eq. (5.6) implementation keeps all five rows regular at the axis:

- `T_{0,k} V_k / X`,
- the radial transport convolution,
- the axial transport convolution,
- radial viscosity, and
- the shifted axial-viscosity term `-Z^[2]_{0,k-1} V_{k-1}/X`.

Differentiating those rows twice in eta requires the existing `RegularFluxThirdMixedJet` plus exactly three additional eta-bearing total-order-four entries:

- `beta_XXetaeta`,
- `beta_Xetaetaeta`, and
- `beta_etaetaetaeta`.

`background_omega_second_parameter_jet.py` exposes these as `RegularFluxFourthMixedJet`.  Its value and first derivative are delegated to the already-landed `omega_over_x_parameter_jet_eq_5_6` after exact projection to the third-mixed jet.  Only the new second derivative is evaluated here, by analytic product/quotient rules.  The nested shifted `Z^[2]` row is differentiated through both `Z` operators without sampling neighboring eta points.

## Verification

`tests/test_background_omega_second_parameter_jet.py` uses independent bivariate polynomial histories of degree two in `X` and degree four in eta, so every newly required fourth-mixed beta entry is nonzero.  For Eq. (5.6) orders `0`, `1`, and `2`, at the regular axis `X=0` and positive `X`, the new bridge must reproduce the already-landed value and first derivative exactly.  Its analytic second derivative is cross-checked against a centered eta difference of the older analytic first-derivative path, with exact polynomial jets rebuilt at the displaced eta points; finite differences are test-only.

Separate sensitivity regressions perturb `beta_XXetaeta`, `beta_Xetaetaeta`, and `beta_etaetaetaeta` independently and verify that each enters the second derivative.  The current-order `beta_XXetaeta` contribution to radial viscosity is also checked separately.  Incomplete histories, weaker/wrong jet types, invalid orders, and non-finite stronger data fail closed.

## Truth boundary

This remains Stage-2 `formal-structure` infrastructure and `paper_exact_velocity_available=false` remains mandatory.

`RegularFluxFourthMixedJet` is an explicit stronger upstream contract.  The current repaired lower-history hierarchy does not yet own these beta fourth-mixed entries; deriving them through Eq. (5.2) requires a stronger U derivative layer.  Caller/test polynomial jets are therefore not paper-exact profile data.  Issue #1 still owns the genuine leading profile.

Accordingly this increment does **not** claim hierarchy-owned `partial_eta^2(Omega/X)`, complete hierarchy-owned `partial_eta^2 actualLowerSource`, `partial_eta^2 f_n`, `partial_eta W_n^(1)`, further Picard convergence, recursively materialized all-order coefficients, the paper's recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity.
