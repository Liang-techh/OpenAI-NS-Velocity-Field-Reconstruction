# Repaired lower-history fourth-mixed provenance

## Scope

This increment connects the already-landed Lemma 5.2 compact fourth-mixed axial repair to the already-landed hierarchy-owned Section 5 lower-history object. It does **not** solve a new coefficient and it does **not** promote the reconstruction to paper-exact status.

## Paper-derived structure

For a positive-order repaired coefficient, the hierarchy now owns one coherent repaired axial source through the fourth-mixed layer required by differentiating Eq. (5.2):

- `U_XXetaeta`,
- `U_Xetaetaeta`,
- `U_etaetaetaeta`.

The fourth-mixed source is constructed only through `Lemma52RepairedFourthMixedJetAdapter`, i.e. the actual compact Lemma 5.2 five-bump repair with analytic `R=sqrt(2X)` derivatives and the supplied finite eta-jets of the repair coefficients. The hierarchy then derives

- the ordinary third-mixed `U` jet by exact projection of that same source,
- the ordinary second jet of `beta=V/X` through the existing analytic Eq. (5.2) adapter, and
- `beta_XXeta`, `beta_Xetaeta`, `beta_etaetaeta` through the existing analytic Eq. (5.2) third-mixed adapter.

No independently supplied beta derivative record is accepted on this path. If a coefficient only owns the older third-mixed U layer, requests for the stronger fourth-mixed U / beta-third-mixed path fail closed.

## Verification

`tests/test_background_repaired_history_fourth_mixed.py` uses an actual `Lemma52MomentRepair` compact geometry and degree-four eta moment data. It cross-checks the hierarchy-owned repaired fourth-mixed U jet against an independently instantiated `Lemma52RepairedFourthMixedJetAdapter`, cross-checks the hierarchy-derived beta third-mixed jet against the adapter's Eq. (5.2) bridge, verifies the ordinary beta-second path is a projection of the same authoritative U source, and verifies absence of fourth-mixed ownership fails closed.

## Truth boundary

This remains `formal-structure` only.

The unrepaired/base `AxialFourthMixedJet`, the moment/patch eta-jets, normalization `C`, and the order-zero leading data remain upstream inputs. In particular, Issue #1 still blocks claiming a fully hierarchy-owned leading fourth-mixed coefficient. Therefore this increment does **not** yet establish hierarchy-owned `partial_eta(Omega/X)` for every strict lower order, hierarchy-owned `partial_eta f_n`, the genuine `k=1` Eq. (5.7) Picard term, Picard convergence, a common analytic strip, the recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity.

`paper_exact_velocity_available=false` remains mandatory.
