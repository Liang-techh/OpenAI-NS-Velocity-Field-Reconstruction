# Hierarchy-owned Eq. (5.6) omega parameter-jet provenance

## Scope

This increment connects two already-landed Stage-2 pieces without solving a new coefficient:

1. `Section5LowerHistoryJetHierarchy` can own coherent fourth-mixed axial `U_j` data and derive the required third-mixed regular-flux `beta_j=V_j/X` jets only through analytic Eq. (5.2).
2. `omega_over_x_parameter_jet_eq_5_6` analytically evaluates `(Omega_k/X, partial_eta(Omega_k/X))` from those beta jets and ordinary axial second jets.

The new bridge accepts only the hierarchy plus `(k,X,eta)`. It does not accept an independent beta derivative table, sampled `Omega/X`, numerical eta derivative, generic cutoff, or fitted coefficient family.

## Paper-derived structure

For one Eq. (5.6) order `k`, the bridge requests coefficient orders `0,...,k` from the contiguous hierarchy. For every such order it:

- obtains the authoritative fourth-mixed `U_j` source,
- derives `beta_j` third-mixed jets through the existing analytic Eq. (5.2) adapter,
- projects the same coherent `U_j` source to the ordinary axial second jet, and
- passes those owned histories to the existing analytic Eq. (5.6) parameter-jet implementation.

If any required coefficient lacks fourth-mixed U ownership, the call fails closed. In particular order zero is not special-cased: until Issue #1 materializes the real leading fourth-mixed profile, this bridge refuses to invent it.

## Verification

`tests/test_background_repaired_history_omega_parameter.py` uses:

- an explicit analytic order-zero fourth-mixed U fixture,
- an actual positive-order `Lemma52RepairedProfileAdapter` built from `Lemma52MomentRepair` compact five-bump geometry and degree-four eta moment data, and
- the hierarchy-owned Eq. (5.2) beta construction.

The returned `Omega_1/X` value is cross-checked against the pre-existing lower-history `actualLowerSource` path. Its analytic eta derivative is then cross-checked independently against a centered eta finite difference of that older value-only path. The regression also checks the axis `X=0` and verifies fail-closed behavior when the leading coefficient does not own fourth-mixed U data.

## Truth boundary

This remains `formal-structure` only and `paper_exact_velocity_available=false` remains mandatory.

The real order-zero fourth-mixed leading profile is still blocked by Issue #1. Positive-order unrepaired/base fourth-mixed U jets, moment/patch eta jets, and normalization data remain upstream inputs. This increment therefore does **not** yet establish a paper-owned complete strict-lower omega-parameter history, hierarchy-owned `partial_eta f_n`, the genuine `k=1` Eq. (5.7) Picard application, Picard convergence, recursive coefficient materialization, the infinite recursive cutoff/local-finiteness argument, Proposition 5.3 all-jets residual decay, or paper-exact velocity.
