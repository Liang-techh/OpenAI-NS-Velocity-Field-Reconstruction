# Eq. (5.7) second-eta Picard primitive provenance

## Scope

This increment adds a reusable exact differentiation primitive for one positive-order Section 5 Picard update

`W^+ = G(A0 W + A1 partial_eta W + f)`.

It is solver infrastructure for Issue #2. It does **not** materialize a new paper coefficient and does not promote the reconstruction truth boundary.

## Analytic identity implemented

Because the singular inverse `G` used in the landed Eq. (5.7) solver commutes with eta differentiation, the primitive assembles

- `W^+ = G(A0 W + A1 W_eta + f)`,
- `W^+_eta = G(A0_eta W + A0 W_eta + A1_eta W_eta + A1 W_etaeta + f_eta)`,
- `W^+_etaeta = G(A0_etaeta W + 2 A0_eta W_eta + A0 W_etaeta + A1_etaeta W_eta + 2 A1_eta W_etaeta + A1 W_etaetaeta + f_etaeta)`.

The last term, `A1 W_etaetaeta`, is an explicit derivative requirement. Therefore the landed hierarchy-owned second-eta matrix and forcing jets are not by themselves sufficient to construct hierarchy-owned `partial_eta^2 W_n^(1)`: the previous Picard iterate must also be owned through its third eta derivative.

## Inputs and truth boundary

The new primitive accepts structural providers for:

- `A0` and `A1` through their second eta derivatives;
- the previous iterate through its third eta derivative;
- forcing through its second eta derivative.

These providers may be analytic test fixtures or future hierarchy-owned adapters. **Caller-supplied providers are not paper-exact merely because they satisfy this interface.** No generic cutoff, finite-order fit, sampled derivative table, or production eta finite difference is introduced.

The repository remains at Stage 2 `formal-structure`, with `full_reconstruction=false` and `paper_exact_velocity_available=false`. Genuine leading-profile ownership remains upstream in Issue #1.

## Verification

`tests/test_background_picard_second_parameter_jet.py` uses independent polynomial matrix/iterate/forcing families with an actual nonzero cubic eta term in the previous iterate. It checks:

1. the new value against the previously landed `picard_map_eq_5_7`;
2. the analytic first and second eta outputs against centered-eta differences of that landed value path (test oracle only);
3. a targeted perturbation of only `W_etaetaeta`, proving that only the second output derivative changes and that the change is exactly `G(A1 delta W_etaetaeta)`;
4. fail-closed validation at `xi=0`, where the singular inverse would otherwise return zero without touching the malformed third-eta input.

No finite difference is used in production code.

## Integration baseline

This increment was branched from main commit `e42cce1c499f7a2e97a303a2fc22c8d06aceb553`, the squash merge of PR #270 (`Add hierarchy-owned second eta matrix jets`).
