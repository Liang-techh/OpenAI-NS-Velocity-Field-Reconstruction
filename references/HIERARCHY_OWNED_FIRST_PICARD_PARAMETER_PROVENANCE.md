# Hierarchy-owned first Picard parameter provenance

## Scope

This increment closes the composition layer between the already landed
hierarchy-owned analytic Eq. (5.7) forcing jet
`(f_n, partial_eta f_n)` and the already landed singular-inverse identity
`partial_eta(G f_n) = G(partial_eta f_n)`.

The implementation is
`src/openai_ns_reconstruction/background_repaired_history_first_picard_parameter.py`.
It accepts only a `Section5LowerHistoryPhiThirdMixedHierarchy`; it does not
accept an independent forcing callback, derivative callback, sampled table,
generic cutoff, or finite-order fitted coefficient.

## Paper connection

Section 5 writes the positive-order system as Eq. (5.7) and Lemma 5.1 uses the
singular diagonal inverse

`(G g)_i(xi,eta) = integral_0^xi (s/xi)^c_i g_i(s,eta) ds`,

with `c = (0,0,2,0,3,1)`.  The kernel is independent of `eta`, so the landed
analytic source jet gives

`partial_eta(G f_n) = G(partial_eta f_n)`.

PR #179 composes the hierarchy-owned forcing provider with the existing
`first_picard_parameter_jet_eq_5_7` primitive.  The resulting value is the
actual first Picard-series term `W_n^(0)=G f_n`, and its parameter component is
the corresponding analytic `partial_eta W_n^(0)`.

This is a real positive-order solver step, but only the `k=0` Picard term.  It
does not yet apply the nontrivial operator
`G(A0 W^(0) + A1 partial_eta W^(0) + f_n)`.

## Fail-closed boundary

The singular inverse returns exactly zero at `xi=0`.  A naive composition could
therefore hide missing Issue-#1 leading mixed data by never evaluating the
forcing at the axis.  The hierarchy bridge explicitly preflights the
hierarchy-owned forcing jet at `xi=0` before applying `G`, so an incomplete
strong leading history still fails closed even for a zero-radius request.

The following remain open:

- the genuine Issue #1 leading mixed profile and its strong derivative data;
- unrepaired positive-order base mixed jets and moment/patch eta-jets;
- the first nontrivial `k=1` Eq. (5.7) Picard application using the landed
  paper-derived `A0/A1` fields;
- Picard convergence and recursive positive-order coefficient materialization;
- the paper's recursively chosen cutoff-scale schedule;
- Proposition 5.3 all-jets truncation/residual decay;
- paper-exact background velocity.

Accordingly Stage 2 remains `formal-structure` and
`paper_exact_velocity_available=false`.

## Validation

`tests/test_background_repaired_history_source_parameter.py` reuses the same
strong order-zero analytic fixture and actual positive-order
`Lemma52MomentRepair` compact five-bump repair used by the preceding hierarchy
bridges.  At recursive order two it checks both the axis and an active compact
repair point.

The value component is cross-checked against the previously landed
`first_picard_term_eq_5_7`, driven by the independently assembled landed
positive-axis forcing value path.  The parameter component is cross-checked
against a centered finite difference in `eta` of that pre-existing first-Picard
value path.  Finite differences appear only in the test oracle; production
applies `G` directly to the analytic hierarchy-owned forcing derivative.

A separate regression removes the leading fourth-mixed U provider and requests
the new bridge at `xi=0`; the request must fail rather than silently returning
`G(0)=0`.  These are floating-point implementation checks, not a proof of
Picard convergence, recursive cutoff admissibility, or all-order residual
decay.
