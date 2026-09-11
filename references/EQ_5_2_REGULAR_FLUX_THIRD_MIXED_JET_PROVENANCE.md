# Section 5 Eq. (5.2) regular-flux third-mixed jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed location: Section 5, Eq. (5.2), used by the positive-axis system in Eqs. (5.3)-(5.7) and Lemma 5.1.  The existing executable second-jet map is `background_regular_flux_second_jets.py`; `background_omega_parameter_jet.py` records the downstream need for three third-mixed regular-flux derivatives.

## Executable mapping

`src/openai_ns_reconstruction/background_regular_flux_third_mixed_jets.py` extends the landed Eq. (5.2) regular-flux map from the ordinary second jet of `beta_n=V_n/X` to the three additional derivatives required by the analytic Eq. (5.6) eta derivative:

- `beta_XXeta`,
- `beta_Xetaeta`, and
- `beta_etaetaeta`.

Eq. (5.2) has numerator

`N = 2 eta U_n - 2(D+lambda_n) eta A[U_n] - (1-eta^2) A[partial_eta U_n]`

and denominator `L=1-2h eta^2`.  Differentiating the previously landed beta second jet once more exposes exactly three new axial-profile derivatives inside the averaged `A[partial_eta U]` term: `U_XXetaeta`, `U_Xetaetaeta`, and `U_etaetaetaeta`.  `AxialFourthMixedJet` makes those requirements explicit instead of manufacturing them numerically.  Radial averages use the exact chain-rule weights `s^2`, `s`, and `1`; the axis `X=0` is therefore evaluated directly, never by sampled division through `V/X`.

The adapter delegates all value/first/second beta fields to the existing `regular_flux_second_jet_eq_5_2` implementation and only adds the three new derivative formulas.  Production contains no neighboring-eta finite difference.

## Independent regression

`tests/test_background_regular_flux_third_mixed_jets.py` uses a polynomial axial family with independent `X^2 eta^2`, `X eta^3`, and `eta^4` terms, so each newly required fourth-mixed U derivative is nonzero.  For orders 0, 1, and 2 and for both `X=0` and positive X, it compares the analytic new beta fields against a centered eta finite difference of the pre-existing Eq. (5.2) second-jet implementation:

- `d_eta(beta_XX)` against `beta_XXeta`,
- `d_eta(beta_Xeta)` against `beta_Xetaeta`, and
- `d_eta(beta_etaeta)` against `beta_etaetaeta`.

A separate regression perturbs only each of the three new fourth-mixed U fields and verifies the corresponding beta derivative changes.  Missing fourth-mixed data, invalid order/domain data, and nonfinite jet entries fail closed.

## Truth boundary

This increment is **formal-structure / solver infrastructure**.  `AxialFourthMixedJet` is still caller/upstream supplied.  The actual Lemma 5.2 compact repair currently materializes U only through the third-mixed jet needed for the beta second jet, and `Section5LowerHistoryJetHierarchy` therefore still owns beta only through second order.

Consequently PR #133 does **not** establish a hierarchy-owned `partial_eta(Omega/X)`, `partial_eta f_n`, the genuine `k=1` Eq. (5.7) Picard term, converged positive-order coefficients, the infinite recursive cutoff schedule, Proposition 5.3 all-jets residual decay, or a paper-exact background.  `paper_exact_velocity_available` remains false.

The strict next bridge is to lift the actual Lemma 5.2 repaired U coefficient to the three fourth-mixed fields above, then let the hierarchy derive this `RegularFluxThirdMixedJet` and feed the already-landed `omega_over_x_parameter_jet_eq_5_6` path.
