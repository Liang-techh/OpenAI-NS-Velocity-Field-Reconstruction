# Section 5 retained recurrence cancellation provenance

## Scope

`background_truncation_residual.py` already implements the exact finite
SlowExpansionResidual recurrence/truncation split and a first-omitted-slow-order
majorant.  It intentionally does **not** infer `R_n=0` from small floating
values.  This increment adds the missing fail-closed admission gate in
`background_recurrence_cancellation.py`.

For every retained order `n=0,...,N`, the hierarchy must supply an exact formal
identity for

`L_n + sum_{i+j=n} K_ij - previous(A)_n = 0`.

Formal atoms carry coefficients in `fractions.Fraction`; floats are rejected.
The three sides are canonically combined and the constructor succeeds only if
the resulting formal expression is exactly empty.  A residual such as
`1/10^30` therefore fails rather than passing a tolerance.  Every identity also
binds one hierarchy id, source revision, coefficient-state id, source/theorem
name, and a nonfuture dependency chain containing its own coefficient order.
The finite certificate requires the complete contiguous prefix `0,...,N` on
one state/revision.

Only after that certificate exists may
`certify_full_residual_tail_majorant(...)` expose the landed omitted-tail
majorant as the conditional full finite residual majorant.  The first omitted
slow order is then exactly `N+1`; the numerical pair/shifted prefactor remains
explicit.

## Truth boundary

This is a finite-prefix exact algebra/provenance gate, not the missing paper
hierarchy and not Proposition 5.3.  The current repository still lacks an
actual hierarchy provider that emits these formal recurrence identities for
arbitrary order, uniform hierarchy-derived `C[j,m]` bounds, and the infinite
recursive cutoff/all-jets-flat argument.  Consequently this layer is hard-coded
`paper_exact=False`; `full_reconstruction=false` and
`paper_exact_velocity_available=false` remain mandatory.

No sampled residual, numerical near-zero decision, tolerance, fitted
coefficient, or generic cutoff enters the exact cancellation decision.

## Verification

`tests/test_background_recurrence_cancellation.py` checks exact rational
cancellation, rejection of a tiny but nonzero rational residual, rejection of
float coefficients, own-order/nonfuture dependency rules, contiguous
same-state prefix ownership, and the gated first-omitted-order majorant.  The
test algebra is a fixture only and is not paper coefficient data.

## Upstream/downstream handoff

Agent 2 remains responsible for the genuine hierarchy-owned coefficient/jets
that define the retained recurrence rows.  Agent 7 must next bind those actual
rows into this exact certificate, combine them with hierarchy-derived
`C[j,m]`, and extend the finite-prefix first-omitted-order statement to the
infinite SlowBorel/DiagonalScale all-jets-flat/super-algebraic closure.
