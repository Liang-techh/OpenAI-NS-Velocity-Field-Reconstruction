# Hierarchy-owned second Picard eta2 provenance

## Scope

This increment closes one finite PositiveAxis/Picard derivative seam in Section 5.
The landed hierarchy path already owns the first Picard-series term
`W_n^(0)=G f_n` through `partial_eta^3` and owns the displayed Eq. (5.7)
matrices `A0/A1` through `partial_eta^2`.  The new bridge constructs the next
Picard-series term

`W_n^(1) = K W_n^(0) = G(A0 W_n^(0) + A1 partial_eta W_n^(0))`

through `partial_eta^2` without accepting caller-maintained matrix, forcing, or
source derivative tables.

## Paper location and identity

Lemma 5.1 / Eq. (5.7) gives

`partial_xi W_n + xi^-1 diag(0,0,2,0,3,1) W_n
 = A0 W_n + A1 partial_eta W_n + f_n`

and the Picard series

`W_n = sum_{k>=0} K^k G f_n`,
`K W = G(A0 W + A1 partial_eta W)`.

Writing `W0=W_n^(0)` and `R=A0 W0 + A1 W0_eta`, exact differentiation gives

- `R_eta = A0_eta W0 + A0 W0_eta + A1_eta W0_eta + A1 W0_etaeta`,
- `R_etaeta = A0_etaeta W0 + 2 A0_eta W0_eta + A0 W0_etaeta
  + A1_etaeta W0_eta + 2 A1_eta W0_etaeta + A1 W0_etaetaeta`.

The singular inverse `G` has kernel `(s/xi)^c_i` and radial endpoints that are
independent of eta, so the finite analytic eta derivatives commute with `G`.
The constructor therefore returns `G R`, `G R_eta`, and `G R_etaeta`.
There is deliberately no extra `f_n` term in this step because `G f_n` is the
separate zeroth Picard-series term `W_n^(0)`.

## Ownership chain

The public constructor accepts only
`Section5LowerHistorySixthMixedHierarchy`.  At every radial quadrature point it
regenerates:

- `W_n^(0)` and eta derivatives 0..3 from the hierarchy-owned first-Picard
  bridge, whose forcing comes from the same strict-lower repaired hierarchy;
- `A0/A1` and eta derivatives 0..2 from the hierarchy-owned PositiveAxis
  matrix bridge, whose base rows come from coefficient-order-zero repaired
  `phi_0`, repaired `U_0`, and analytic Eq. (5.2) `beta_0=V_0/X`.

No caller `SourceJet`, repair jet, forcing derivative table, matrix derivative
table, independent normalization constant, sampled `V/X`, generic cutoff, or
fitted coefficient is admitted.  No new seventh-mixed U tier is introduced.

## Verification

Focused regression uses a coherent analytic strong-hierarchy fixture and checks:

- the value row is exactly `G(A0 W0 + A1 W0_eta)`;
- the second eta row equals an independently assembled exact second-Leibniz
  right-hand side sent through the same manuscript singular inverse;
- all three rows vanish exactly at the axis and are read-only;
- wrong hierarchy types fail closed.

The regression uses no finite-difference oracle.  Its analytic fixture is test
data only and is not promoted to paper coefficient data.

## Truth boundary

Status remains Stage-2 `formal-structure`.

This increment constructs only the finite Picard-series term `W_n^(1)` through
eta order two.  It does **not** establish convergence of the Picard series, a
repaired positive-order `(phi_n,U_n,Pi_n,V_n)`, a total arbitrary-order
coefficient provider, theorem-backed `C[j,m]` rows, the recursive infinite
cutoff/all-jets-flat closure, Proposition 5.3, full reconstruction, or
paper-exact velocity.

Uniform/stagewise bounds, recursive cutoff scheduling, and finite-to-all-order
residual convergence remain the separate Agent-7 workstream.  The next
Agent-2 coefficient seam is to continue the hierarchy-owned Picard recursion
far enough to materialize an actual positive-order `W_n`, then extract and
repair its `phi_n/U_n/Pi_n/V_n` coefficient data under the same hierarchy
ownership rather than returning to caller-supplied jets.
