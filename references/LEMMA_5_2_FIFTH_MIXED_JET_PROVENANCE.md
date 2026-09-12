# Section 5 Lemma 5.2 repaired fifth-mixed U-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed locations: Section 5, Lemma 5.2 and Eqs. (5.14)-(5.16) for the compact five-bump moment repair, Eq. (5.2) for `beta_n=V_n/X`, and Eqs. (5.3)-(5.7) for the positive-axis recursion.

## Executable mapping

`src/openai_ns_reconstruction/background_moment_repair_fifth_mixed_jets.py` lifts the already-landed compact Lemma 5.2 repaired U jet by exactly the three total-order-five derivatives required by `background_regular_flux_fourth_mixed_jets.py`:

- `U_XXetaetaeta`,
- `U_Xetaetaetaeta`, and
- `U_etaetaetaetaeta`.

For `Delta U_n(X,eta)=sum_j alpha_j(eta)b_j(sqrt(2X))`, these are evaluated analytically as `sum_j alpha_j''' d_X^2 b_j`, `sum_j alpha_j'''' d_X b_j`, and `sum_j alpha_j''''' b_j`. The same compact five bumps and the same Lemma 5.2 coefficient solve are used; no pointwise toy repair is introduced.

The first twelve fields are delegated to `Lemma52RepairedFourthMixedJetAdapter` using the exact fourth-order projection of the same caller-supplied fifth-mixed base provider. This prevents the stronger layer from silently disagreeing with the already-landed repaired fourth-mixed history. Off compact support the unrepaired fifth-mixed jet is returned exactly and unavailable fifth eta rows are not requested.

`beta_fourth_mixed_jet(...)` feeds this repaired U jet directly into the landed analytic Eq. (5.2) fourth-mixed regular-flux map, retaining its positive-order `lambda_n=2nh` correction rather than reimplementing Eq. (5.2).

## Independent regression

`tests/test_background_moment_repair_fifth_mixed_jets.py` uses the actual `Lemma52MomentRepair` five-bump geometry with degree-five analytic moment functions. The three new repaired U entries are checked against centered eta differences of the previously landed repaired fourth-mixed adapter, while the complete fourth-order projection must agree exactly.

A second end-to-end check feeds the repaired fifth-mixed U jet through the landed Eq. (5.2) fourth-mixed beta adapter and compares its three new beta derivatives against centered eta differences of the previous repaired beta third-mixed path. Finite differences are test-only. Regression also verifies compact off-support exactness, active-support fail-closed behavior when fifth eta rows are absent, and rejection of a non-`AxialFifthMixedJet` base provider.

## Truth boundary

This increment remains **formal-structure / solver infrastructure**. The unrepaired/base `AxialFifthMixedJet` and the moment/patch eta jets remain upstream inputs. They are not inferred from sampled data and are not promoted to paper-exact status. In particular, Issue #1 still owns the genuine leading profile data.

This bridge supplies the real compact-repair contribution needed for a hierarchy-owned `RegularFluxFourthMixedJet`, but the Section 5 lower-history hierarchy still has to own this stronger provider and route it into the landed second-eta Eq. (5.6) `Omega/X` path. Hierarchy-owned `partial_eta^2 actualLowerSource`, `partial_eta^2 f_n`, `partial_eta W_n^(1)`, Picard convergence, recursive coefficient materialization, the recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, and paper-exact velocity remain open.

`paper_exact_velocity_available` remains false.
