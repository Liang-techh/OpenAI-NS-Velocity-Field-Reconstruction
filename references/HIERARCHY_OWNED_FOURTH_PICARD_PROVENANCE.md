# Hierarchy-owned fourth Picard value provenance

## Scope

This increment advances the Stage-2 Section-5 PositiveAxis recursion by one
bounded Picard-series term. For positive coefficient order `n`, Lemma 5.1 and
Eq. (5.7) write

`W_n = sum_{k>=0} K^k G f_n`

with

`K W = G(A0 W + A1 partial_eta W)`.

The landed stacked hierarchy path owns the finite derivative triangle

`W_n^(0): eta^3 -> W_n^(1): eta^2 -> W_n^(2): eta`.

This increment consumes the final available derivative row and constructs

`W_n^(3)=K W_n^(2)`

at value level.

## Source binding

Pinned formal source:

- `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/PositiveAxisSystem.lean`
- paper Lemma 5.1 / Eq. (5.7)

The state vector is the pinned six-vector
`(phi_n, U_n, K_n, Pi_n, partial_xi phi_n, partial_xi U_n)`, and `G` is the
already-landed componentwise singular inverse from Eq. (5.7).

## Ownership chain

Production accepts only `Section5LowerHistorySixthMixedHierarchy`. At every
radial quadrature point it regenerates:

1. hierarchy-owned `W_n^(2)` value and first eta derivative from the preceding
   Picard bridge;
2. hierarchy-owned `A0/A1` values from coefficient-order-zero strong profile
   data;
3. the fourth Picard right-hand side
   `A0 W_n^(2) + A1 partial_eta W_n^(2)`;
4. the output value by applying the same singular inverse `G`.

There is no extra forcing term in `W_n^(3)`: `G f_n` is the separate
`W_n^(0)` term of the Lemma 5.1 series.

## Forbidden substitutions

The production path accepts no caller-supplied `SourceJet`, repaired jet,
forcing derivative table, matrix derivative table, or independent
normalization data. It uses no finite differences, fitted derivatives, generic
cutoff, sampled `V/X` division, or finite residual tolerance.

## Verification

The focused regression independently regenerates hierarchy-owned `W_n^(2)` and
`A0/A1` at the singular-inverse quadrature points, forms the displayed Eq.
(5.7) right-hand side, and reapplies the landed `G`. It also checks exact axis
zero behavior, read-only output, finite shape, and fail-closed hierarchy
typing. The analytic fixture is regression data only and is not paper
coefficient data.

## Workstream split

Agent 7 continues to own stagewise/exact-majorant `C[j,m]`, recursive
SlowBorel/DiagonalScale scheduling, common support/uncut-radius evidence,
retained-recurrence evidence, target-order arithmetic, and finite-to-all-order
convergence gates. This increment does not implement or infer those results.

## Truth boundary

Status remains `formal-structure`; `full_reconstruction=false` and
`paper_exact_velocity_available=false`.

Completing this finite eta-derivative triangle does **not** prove convergence
of the Lemma 5.1 Picard series and does not materialize a repaired positive-
order `(phi_n,U_n,Pi_n,V_n)`. It does not provide a total arbitrary-order
coefficient provider, recursive-cutoff/all-jets-flat closure, Proposition 5.3,
full reconstruction, or paper-exact velocity. The next Stage-2 step must move
toward a genuine Picard-tail/contraction or coefficient-extraction bridge
rather than treating this finite prefix as a solved coefficient.
