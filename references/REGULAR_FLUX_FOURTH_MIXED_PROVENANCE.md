# Section 5 Eq. (5.2) fourth-mixed regular-flux provenance

## Scope

The landed `background_omega_second_parameter_jet.py` evaluates the analytic second eta derivative of the regular Eq. (5.6) quotient `Omega_k/X`.  Its input contract is a fourth-mixed jet of `beta_k = V_k/X`.  Before this increment that stronger beta jet could only be supplied by a caller.  This increment adds the analytic Eq. (5.2) derivative map needed to construct that beta layer from a stronger axial `U_k` jet.

## Paper-derived structure

Eq. (5.2) is used in the already-landed regular form

`L beta_n = 2 eta U_n - 2 c_n eta A[U_n] - (1-eta^2) A[partial_eta U_n]`,

with `c_n = 1/2-h+lambda_n` and the positive-order `lambda_n = 2nh` correction retained by the existing adapters.  `background_regular_flux_fourth_mixed_jets.py` does not replace those lower-order adapters.  It delegates the complete third-mixed beta jet to `regular_flux_third_mixed_jet_eq_5_2` and differentiates only the three additional entries required by the landed Eq. (5.6) second-eta path:

- `beta_XXetaeta`,
- `beta_Xetaetaeta`, and
- `beta_etaetaetaeta`.

Because Eq. (5.2) contains `A[partial_eta U]`, these require exactly three additional eta-bearing fifth-mixed axial entries:

- `U_XXetaetaeta`,
- `U_Xetaetaetaeta`, and
- `U_etaetaetaetaeta`.

They are represented by `AxialFifthMixedJet`.  Radial averages are differentiated analytically with chain-rule weights `s^2`, `s`, and `1`.  The axis `X=0` is therefore handled directly; production never forms a sampled `V/X` quotient and never finite-differences in eta.

## Verification

`tests/test_background_regular_flux_fourth_mixed_jets.py` uses an independent bivariate polynomial axial profile of degree two in `X` and degree five in eta, chosen so all three newly required fifth-mixed U entries are nonzero.  For coefficient orders `0`, `1`, and `2`, at `X=0` and positive `X`, the new analytic fourth-mixed beta entries are cross-checked against centered eta differences of the already-landed analytic third-mixed Eq. (5.2) path.  Finite differences are test-only.

The test also requires the lower projection `result.third()` to equal the existing third-mixed result exactly, perturbs each new fifth-mixed U entry independently to verify structural activity, and checks fail-closed behavior for weaker fourth-mixed providers, invalid domains, invalid orders, and non-finite stronger data.

## Truth boundary

This remains Stage-2 `formal-structure` infrastructure and `paper_exact_velocity_available=false` remains mandatory.

`AxialFifthMixedJet` is an explicit stronger upstream contract.  The Lemma 5.2 repaired hierarchy currently owns repaired U only through the fourth-mixed layer, so this increment does **not** yet make the Eq. (5.6) second-eta input hierarchy-owned.  Caller/test polynomial fifth-mixed jets are not paper-exact profile data, and Issue #1 still owns the genuine leading profile.

Accordingly this increment does **not** claim a paper-exact positive-order coefficient, hierarchy-owned `partial_eta^2(Omega/X)`, complete hierarchy-owned `partial_eta^2 actualLowerSource` or `partial_eta^2 f_n`, `partial_eta W_n^(1)`, Picard convergence, recursive coefficient materialization, the recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity.