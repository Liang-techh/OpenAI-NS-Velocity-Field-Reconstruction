# Section 5 Eq. (5.6) Omega/X parameter-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: Section 5, especially Eqs. (5.3)-(5.7) and Lemma 5.1.  The executable value formula is already pinned in `background_recurrence.py` from the repository's `PositiveAxisSystem` formalization.

## Executable mapping

`src/openai_ns_reconstruction/background_omega_parameter_jet.py` differentiates the landed regular Eq. (5.6) expression for `Omega_k/X` once in eta.  It covers all five rows in that expression: the scaling operator, radial transport convolution, axial transport convolution, radial viscosity, and the shifted nested axial-viscosity operator.

The production derivative is analytic.  In particular, the nested shifted axial-viscosity row differentiates both regular `Z` operators with quotient/product rules and requires the explicit third mixed regular-flux derivatives `beta_XXeta`, `beta_Xetaeta`, and `beta_etaetaeta`.  No neighboring-eta finite difference is used in production.

`tests/test_background_omega_parameter_jet.py` uses polynomial beta/U families containing `X^2 eta`, `X eta^2`, and `eta^3` terms.  It compares the analytic derivative with an independent centered finite difference of the pre-existing `omega_over_x_eq_5_6` value path for orders 0, 1, and 2, including the axis `X=0`.  A separate regression perturbs only the three new third-mixed beta entries to ensure those inputs are not silently ignored.

## Truth boundary

This increment is **formal-structure / solver infrastructure**.  The current `Section5LowerHistoryJetHierarchy` owns beta only through second order, derived analytically from Eq. (5.2).  Constructing the third-mixed beta jet required here needs one additional derivative layer of the axial coefficient U (and corresponding Lemma 5.2 repair jets).  Until that data is materialized, this module must not be described as a hierarchy-owned `partial_eta f_n` or as a paper-exact positive-order coefficient.

The strict next step is to extend the Eq. (5.2) regular-flux adapter to the three third-mixed beta derivatives with explicit higher-U-jet requirements, then connect this module to `actualLowerSource` and the existing `first_picard_parameter_jet_eq_5_7`.  Only after that can the repository form the genuine `k=1` Picard term from hierarchy-owned data.

No full Picard convergence certificate, recursively materialized all-order coefficient sequence, common analytic strip, infinite cutoff schedule, Proposition 5.3 all-jets residual decay, or final cutoff-summed background is established here.  `paper_exact_velocity_available` remains false.
