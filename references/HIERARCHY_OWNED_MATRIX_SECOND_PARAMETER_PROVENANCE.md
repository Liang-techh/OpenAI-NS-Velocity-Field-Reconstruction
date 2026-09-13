# Hierarchy-owned Eq. (5.7) matrix second-eta jet provenance

## Scope

This increment is a Stage-2 Section-5 `formal-structure` step based on `main`
commit `d875ab912b137a7d909b087ade6faaa12d69ec89`. It adds the analytic second
`eta` derivative of the displayed Eq. (5.7) matrices needed before the
hierarchy-owned Picard chain can be differentiated beyond the landed `k=1`
parameter jet.

The existing bridge already owns `(A0, partial_eta A0)` and
`(A1, partial_eta A1)`. The new path keeps those value/first-derivative outputs
authoritative and computes only `partial_eta^2 A0` and `partial_eta^2 A1`.

## Construction

`background_repaired_history_matrix_second_parameter.py` evaluates the exact
PositiveAxis matrix formulas with an internal second-order scalar eta jet. Its
inputs are all hierarchy-owned:

- the order-zero third-mixed `phi_0` jet;
- the order-zero third-mixed `U_0` jet; and
- the order-zero second jet of `beta_0=V_0/X`, still derived only through the
  landed analytic Eq. (5.2) adapter.

The scalar jet applies exact product/quotient rules to the displayed formulas.
No production finite difference, sampled matrix table, generic cutoff, fitted
coefficient family, or caller-supplied matrix derivative is accepted.

## Verification

`tests/test_background_repaired_history_matrix_second_parameter.py` reuses a
strong analytic leading-history fixture. At the axis and away from it, it
requires the new value and first derivative matrices to be exactly identical to
the landed matrix-parameter bridge. The new second derivatives are independently
cross-checked against centered eta differentiation of that already-landed
analytic first-derivative bridge. The finite difference is test-only.

The same regression checks fail-closed behavior for an unowned recursive order
and for a hierarchy that does not provide the required strong type.

## Truth boundary

This remains `formal-structure` only:

- genuine Issue-#1 leading profile data remain upstream;
- caller/test analytic jets are not manuscript data;
- second matrix eta derivatives do not establish a converged Picard series;
- no finalized recursive coefficient family is materialized;
- no paper recursive cutoff-scale schedule is claimed;
- no Proposition 5.3 arbitrary-order truncation/residual-decay claim is made;
- `full_reconstruction=false`; and
- `paper_exact_velocity_available=false`.
