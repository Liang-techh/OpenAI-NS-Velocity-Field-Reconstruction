# Hierarchy-owned third Picard first-parameter provenance

## Scope

This increment advances the Stage-2 Section-5 PositiveAxis recursion by one
bounded Picard-series term. For positive coefficient order `n`, Lemma 5.1 and
Eq. (5.7) write

`W_n = sum_{k>=0} K^k G f_n`

with

`K W = G(A0 W + A1 partial_eta W)`.

The landed hierarchy path already owns `W_n^(1)=K W_n^(0)` through
`partial_eta^2`. This increment constructs

`W_n^(2)=K W_n^(1)`

through `partial_eta`.

## Source binding

Pinned formal source:

- `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/PositiveAxisSystem.lean`
- paper Lemma 5.1 / Eq. (5.7)

The state vector is the pinned six-vector
`(phi_n, U_n, K_n, Pi_n, partial_xi phi_n, partial_xi U_n)`, and the singular
inverse is the already-landed componentwise Eq. (5.7) operator `G`.

## Ownership chain

Production accepts only `Section5LowerHistorySixthMixedHierarchy`. At every
radial quadrature point it regenerates:

1. hierarchy-owned `W_n^(1)` value / eta / eta^2 rows;
2. hierarchy-owned `A0/A1` value / eta rows from coefficient-order-zero strong
   profile data;
3. the third Picard right-hand side
   `A0 W1 + A1 partial_eta W1`;
4. its exact first eta derivative
   `A0_eta W1 + A0 W1_eta + A1_eta W1_eta + A1 W1_etaeta`;
5. the two output rows by applying the same singular inverse `G`.

There is no extra forcing term in `W_n^(2)`: `G f_n` is the separate
`W_n^(0)` term.

## Forbidden substitutions

The production path does not accept caller-supplied `SourceJet`, repaired jet,
forcing derivative table, matrix derivative table, or independent
normalization data. It uses no finite differences, fitted derivatives,
generic cutoff, sampled `V/X` division, or finite residual tolerance.

## Verification

The focused regression independently reassembles the value and first-eta
right-hand sides from the hierarchy-owned `W_n^(1)` and matrix jets, applies
the landed singular inverse, checks exact axis-zero behavior/read-only output,
and verifies fail-closed hierarchy typing. The analytic fixture is regression
data only and is not paper coefficient data.

## Truth boundary

Status remains `formal-structure`.

This increment does **not** prove convergence of the Picard series, does not
materialize a repaired positive-order `(phi_n,U_n,Pi_n,V_n)`, does not provide
a total arbitrary-order coefficient provider, and does not establish the
recursive cutoff/all-jets-flat argument, Proposition 5.3, full reconstruction,
or paper-exact velocity. Agent 7 continues to own stagewise/exact-majorant
`C[j,m]`, recursive SlowBorel/DiagonalScale scheduling, common support,
retained-recurrence evidence, target-order arithmetic, and convergence gates.
