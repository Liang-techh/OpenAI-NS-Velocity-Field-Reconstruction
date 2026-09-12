# Hierarchy-owned nontrivial k=1 Picard provenance

## Scope

PR #184 closes the next Eq. (5.7) composition layer after PR #179.  The new
implementation is
`src/openai_ns_reconstruction/background_repaired_history_k1_picard.py` and
accepts only a `Section5LowerHistoryPhiThirdMixedHierarchy`.

It evaluates

`W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n)`

with every input obtained from already-landed hierarchy-owned infrastructure:

- `A0` and `A1` come only from `hierarchy.positive_axis_fields(order)`, which
  constructs the displayed PositiveAxis matrices from the owned base jets;
- `(W_n^(0), partial_eta W_n^(0))` comes only from PR #179's strong-history
  first-Picard parameter bridge;
- `f_n` comes only from the hierarchy-owned analytic forcing bridge landed in
  PR #174; and
- the radial inverse is the existing Eq. (5.7) singular inverse used by the
  canonical `picard_map_eq_5_7` primitive.

The interface accepts no caller-supplied matrices, forcing, previous iterate,
previous eta derivative, sampled derivative table, generic cutoff, or fitted
coefficient.

## Paper connection

Section 5 / Lemma 5.1 rewrites the positive-order coefficient system as

`partial_xi W_n + xi^-1 diag(0,0,2,0,3,1) W_n = A0 W_n + A1 partial_eta W_n + f_n`

and applies the singular inverse `G` in a Picard construction.  PR #179 landed
the hierarchy-owned `k=0` value/parameter jet `W_n^(0)=G f_n`.  PR #184 now
performs the first nontrivial application of the full displayed right-hand side
operator to that iterate.

This value is a genuine hierarchy-owned Picard iterate on the currently owned
strict-lower history.  It is not the isolated series term `K G f_n`, not the
converged fixed point, and not a finalized recursive coefficient.

## Fail-closed boundary

The low-level singular inverse returns exactly zero at `xi=0`, so a direct
wrapper could otherwise skip evaluation of the strong history.  The PR #184
bridge preflights the PR #179 first-Picard parameter jet, the hierarchy-owned
forcing, and the structured `A0/A1` matrices at the axis before calling the
Picard map.  Missing Issue-#1 fourth-mixed leading data therefore still fail
closed even on a zero-radius request.

The order-zero mixed profile used in tests is explicitly an analytic fixture.
The genuine leading profile remains an Issue #1 upstream dependency.  Likewise,
unrepaired positive-order base mixed jets, moment/patch eta-jets, the recursive
cutoff-scale schedule, Picard convergence, coefficient materialization, and
Proposition 5.3 truncation/residual decay are not supplied by this increment.

Accordingly Stage 2 remains `formal-structure`, `full_reconstruction=false`,
and `paper_exact_velocity_available=false`.

## Validation

`tests/test_background_repaired_history_source_parameter.py` reuses the existing
real `Lemma52MomentRepair` compact five-bump repaired order-1 coefficient and
runs the new bridge at recursive order 2.  At an off-axis point in the compact
repair regime, the returned value is cross-checked against an independently
assembled right-hand side

`A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n`

passed directly to the pre-existing singular inverse.  The test additionally
checks that this nontrivial iterate is not numerically identical to the landed
`k=0` iterate for the fixture.

Separate regressions verify exact zero on the axis only after successful strong
preflight, reject a missing leading fourth-mixed U provider even at `xi=0`, and
reject non-strong hierarchy objects.  These are floating-point implementation
checks, not a proof of Picard convergence or paper-exact reconstruction.
