# Hierarchy-owned PositiveAxis matrix eta-jet provenance

## Scope

This increment closes one finite derivative-ownership seam in the Section 5
Eq. (5.7) PositiveAxis system.  The landed hierarchy bridge already constructs
the value matrices `A0` and `A1` from coefficient order zero.  The new bridge
constructs those same matrices through `partial_eta^2` without accepting a
caller-maintained matrix derivative table.

## Paper location and identity

Section 5, Eq. (5.7), has

`partial_xi W_n + xi^-1 diag(0,0,2,0,3,1) W_n
 = A0 W_n + A1 partial_eta W_n + f_n`.

The displayed `A0/A1` entries use only the order-zero leading quantities
`phi_0`, `U_0`, `beta_0=V_0/X`, their displayed X/eta derivatives, and the
explicit geometry factors `1-eta^2` and `ell=1-2h eta^2`.

`background_repaired_history_positive_axis_matrix_jets.py` differentiates
exactly those landed formulas through eta order two.  Scalar products use the
ordinary first/second Leibniz rule.  Quotients use the exact recurrence from
`b q = a`:

- `q = a/b`,
- `q' = (a' - b' q)/b`,
- `q'' = (a'' - 2 b' q' - b'' q)/b`.

No finite difference or fitted derivative is used in production.

## Ownership chain

The public constructor accepts only
`Section5LowerHistorySixthMixedHierarchy`.  At `X=xi^2` it obtains:

- `phi_0` from the hierarchy-owned repaired fifth-mixed phi jet,
- `U_0` from the hierarchy-owned repaired sixth-mixed axial jet,
- `beta_0` from the same `U_0` through the analytic Eq. (5.2) adapter.

The existing hierarchy-owned matrix constructor remains authoritative for the
value rows.  The eta derivatives consume only exact projections of the same
order-zero strong hierarchy.  No independent `PositiveAxisBaseJet`, beta table,
matrix jet, `SourceJet`, repair jet, sampled `V/X`, generic cutoff, or fitted
coefficient is admitted.

The existing fifth/sixth mixed hierarchy is already stronger than required:
`partial_eta^2 A0/A1` needs at most third eta derivatives of `phi_0/U_0`,
second eta derivatives of their X rows, and second eta derivative of `beta_0`.
No seventh-mixed U layer is introduced.

## Verification

Focused regression uses analytic polynomial strong-jet fixtures.  It verifies
exact value-row delegation and compares the new derivatives with an independent
test-only normalized Taylor-series algebra whose degree-two coefficient is
`f''/2`.  The oracle uses exact finite series multiplication/division rather
than finite differences.  Axis finiteness, read-only outputs, and fail-closed
hierarchy typing are also checked.

## Truth boundary

Status remains Stage-2 `formal-structure`.

This increment does **not** construct `W_n^(1)`, a converged Picard sum,
a repaired positive-order coefficient `(phi_n,U_n,Pi_n,V_n)`, an arbitrary-order
coefficient provider, uniform/stagewise `C[j,m]`, recursive cutoff convergence,
Proposition 5.3, full reconstruction, or paper-exact velocity.  The next
Agent-2 seam is to combine these matrix eta jets with the landed
`W_n^(0)[0..3]` jet to construct hierarchy-owned `partial_eta^2 W_n^(1)`.
Uniform bounds/cutoff/convergence remain the separate Agent-7 workstream.
