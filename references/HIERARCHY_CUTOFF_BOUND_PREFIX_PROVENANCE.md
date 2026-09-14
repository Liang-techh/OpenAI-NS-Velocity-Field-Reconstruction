# Section 5 hierarchy-bound cutoff-prefix provenance

## Scope

The landed `background_cutoff_schedule.py` implements the pinned
`SlowBorelBase.exists_admissibleScales` / `DiagonalScale.exists_diagonal_scales`
finite-prefix scale selection once normalized-template constants `C[j,m]` are
supplied. Its remaining analytic ownership seam is that a plain callable can
supply those numbers without recording which hierarchy coefficients justify
them.

This increment closes only that provenance/ownership seam. It introduces a
dependency-tracked finite-prefix certificate in
`background_hierarchy_cutoff_bounds.py`: every `C[j,m]` must name one hierarchy
identity, one source revision, and a nonempty coefficient dependency chain that
includes coefficient order `j` and contains no future order. A certificate is
accepted only when it contains the exact triangular set

`1 <= j <= J`, `0 <= m <= j+2`.

Only after those checks may the bounds enter the existing SlowBorel cutoff
scheduler.

The implementation is merge-forwarded onto repository `main`
`51e7b91bf5bdd9dfd52ecc74f3e1b8cb61def464`, which includes Agent 2 PR #301's
hierarchy-owned third-eta Omega/X derivative row. This integration changes no
Agent-2 hierarchy jet or coefficient-recurrence implementation; those remain
upstream inputs to this Agent-7 lane.

## Mathematical boundary

This does **not** derive the paper's uniform `C[j,m]` constants from the current
finite hierarchy. Instead it provides the fail-closed contract that a future
analytic hierarchy-bound layer must satisfy. In particular:

- no sampled derivative maximum is accepted;
- no plain untracked float provider is promoted by this bridge;
- every bound must bind its own coefficient order and may refer only to current
  or lower coefficient orders;
- hierarchy identity and source revision must be constant over the whole
  finite prefix;
- missing, duplicate, wrongly indexed, nonpositive, nonfinite, cross-hierarchy,
  or cross-revision entries fail before schedule construction;
- the resulting witness is hard-coded as non-paper-exact.

The existing scheduler remains authoritative for the exact finite-prefix
`C[j,m] a_j^(-h j) <= 2^(-j)` log-space checks, recursive doubling envelope,
and arbitrary-precision integer scale path.

## Verification

`tests/test_background_hierarchy_cutoff_bounds.py` uses a typed analytic-bound
provider fixture, not paper coefficient data. It checks the exact triangular
SlowBorel index set, preservation of certified rows into the scheduler, every
scheduled dyadic-edge margin, own-order dependency binding, fail-closed malformed
or cross-revision inputs, and that `paper_exact` remains false.

Repository CI is authoritative for integration and the complete suite after
this merge-forward.

## Remaining boundary for Agent 2 / Agent 7

Stage 2 remains `formal-structure`, with `full_reconstruction=false` and
`paper_exact_velocity_available=false`.

Agent 2 must still provide, for every requested hierarchy order, genuine
hierarchy-owned analytic derivative/support rows and exact retained recurrence
row decompositions on one coefficient state. Agent 7 still needs those rows to
be controlled by one all-order envelope strong enough for the pinned
DiagonalGrowth quantifiers: for every derivative budget `m`, the effective gain
`gain(j) - m * cutLoss(j)` must tend to `+infinity`. Only then can the current
finite-prefix certificates be promoted to a constructive arbitrary-order
extension, all-jets-flatness, Proposition 5.3, and super-algebraic residual
decay.

A green test for this finite-prefix bridge is therefore not evidence of
all-order convergence or a paper-exact background.
