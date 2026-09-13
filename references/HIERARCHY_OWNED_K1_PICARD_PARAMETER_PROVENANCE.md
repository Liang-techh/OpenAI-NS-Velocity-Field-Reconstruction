# Hierarchy-owned k=1 Picard parameter-jet provenance

## Scope

This increment composes already-landed Section 5 analytic layers to expose the first eta derivative of the first nontrivial Eq. (5.7) Picard iterate,

`(W_n^(1), partial_eta W_n^(1))`.

It does not introduce a new coefficient ansatz, generic cutoff, sampled derivative table, or numerical fitting stage. The value path remains the already-landed hierarchy-owned `k=1` Picard application; only its analytic eta derivative is new.

## Paper/formula binding

For positive order `n`, Eq. (5.7) is written as

`W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n)`.

The singular inverse `G` has an eta-independent kernel. Differentiating the displayed Picard right-hand side at fixed `xi` therefore gives

`partial_eta W_n^(1) = G(`
`    (partial_eta A0) W_n^(0)`
`  + A0 partial_eta W_n^(0)`
`  + (partial_eta A1) partial_eta W_n^(0)`
`  + A1 partial_eta^2 W_n^(0)`
`  + partial_eta f_n` 
`)`.

Every term in this formula is supplied by a previously landed analytic path:

- `(A0, partial_eta A0)` and `(A1, partial_eta A1)` come from the hierarchy-owned Eq. (5.7) matrix parameter jets;
- `(W_n^(0), partial_eta W_n^(0), partial_eta^2 W_n^(0))` is obtained by applying the existing first-Picard second-parameter primitive to the hierarchy-owned `(f_n, partial_eta f_n, partial_eta^2 f_n)` forcing;
- the value `W_n^(1)` delegates unchanged to the existing hierarchy-owned `k=1` value implementation;
- the new parameter component is evaluated by the same canonical Eq. (5.7) singular inverse.

## Ownership and fail-closed boundary

The public bridge accepts only `Section5LowerHistoryPhiFourthMixedHierarchy`. It accepts no independent matrix, forcing, previous iterate, eta derivative, beta/Omega table, generic cutoff, or fitted coefficient family.

The strong hierarchy therefore remains responsible for the repaired Lemma 5.2 fourth-mixed `phi` and fifth-mixed `U` history from which the landed lower-source/forcing derivative chain is derived. The genuine order-zero fourth/fifth-mixed leading profile is still an Issue #1 dependency. At `xi=0`, the bridge explicitly preflights its strong forcing and matrix dependencies before the exact `G(0)=0` shortcut can return, so missing Issue-#1 data remain fail-closed rather than being hidden by the axis value.

## Regression

`tests/test_background_repaired_history_forcing_second_parameter.py` reuses the repository's actual `Lemma52MomentRepair.from_intervals(...)` five-bump compact repair in a recursive order-2 hierarchy. The new `W_n^(1)` value must be exactly identical to the previously landed hierarchy-owned value path. The analytic eta derivative is cross-checked against centered eta differentiation of that independent landed value path at the axis and at a point inside active compact-repair support. Finite differences are test-only oracles and are not used by production.

The regression also verifies that the new public bridge rejects a hierarchy that does not own the required strong fourth-mixed `phi` layer.

## Truth boundary

Status remains `formal-structure`. `paper_exact_velocity_available=false` and `full_reconstruction=false`.

This increment establishes one hierarchy-owned differentiated Picard iterate only. It does **not** establish Picard convergence, a finalized positive-order coefficient, the complete all-order recursive coefficient family, the paper's recursively chosen cutoff-scale schedule, Proposition 5.3 all-jets truncation/residual decay, or paper-exact velocity. Caller/test leading fourth/fifth-mixed data are not substitutes for the genuine Issue #1 leading profile.
