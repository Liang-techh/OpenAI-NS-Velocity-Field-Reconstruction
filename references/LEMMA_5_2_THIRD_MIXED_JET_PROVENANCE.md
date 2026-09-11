# Lemma 5.2 repaired axial third-mixed jet provenance

Status: **formal-structure**. This record does not promote Stage 2 to paper-exact and does not make `paper_exact_velocity_available` true.

## Paper location and role

This adapter sits between the compact moment repair of Lemma 5.2 / Eq. (5.14)-(5.16) and the Eq. (5.2) regular flux reconstruction. The landed `background_regular_flux_second_jets.py` needs the repaired axial coefficient through `U_XXeta`, `U_Xetaeta`, and `U_etaetaeta` in order to form the full second `(X,eta)` jet of `beta_n=V_n/X`.

For the compact correction

`Delta U(X,eta) = sum_j alpha_j(eta) b_j(sqrt(2X))`,

the production adapter uses the exact product rules

- `d_X^2 d_eta Delta U = sum_j alpha'_j d_X^2 b_j`,
- `d_X d_eta^2 Delta U = sum_j alpha''_j d_X b_j`,
- `d_eta^3 Delta U = sum_j alpha'''_j b_j`.

The already-landed analytic `R=sqrt(2X)` bump derivatives are reused, and `Lem​ma52RepairedProfileAdapter.repair_coefficients(..., max_order=3)` supplies ordinary eta derivatives of the repair coefficients. No production finite difference is used.

## Independent regression boundary

`tests/test_background_moment_repair_third_mixed_jets.py` differentiates only independently evaluable repaired `U_n` function values with centered finite-difference stencils to check the three new third-mixed entries. A separate check compares the resulting Eq. (5.2) `beta_n` value against the existing function-level `LeadingProfile.radial_flux_factor` path. Off the compact repair supports the adapter returns the unrepaired jet exactly and does not demand unavailable third eta moment data.

## What remains unresolved

The unrepaired `AxialThirdMixedJet` and the moment/patch eta-jets are still caller-supplied. This increment therefore does **not** prove that those data come from the genuine converged positive-order hierarchy, does not certify a common analytic strip or uniform nonvanishing patch factor, does not auto-materialize every lower-history order, and does not establish the recursive infinite cutoff sum or Proposition 5.3 all-jets-flat residual decay.

The next recursive integration step is to build hierarchy-owned repaired coefficient objects whose second jets and Eq. (5.2) `beta_n` second jets can be inserted into `actualLowerSource` without caller-assembled derivative records.
