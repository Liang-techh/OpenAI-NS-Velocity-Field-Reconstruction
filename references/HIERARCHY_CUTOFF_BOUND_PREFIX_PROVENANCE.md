# Section 5 hierarchy-bound cutoff-prefix provenance

## Scope

The landed `background_cutoff_schedule.py` implements the pinned
`SlowBorelBase.exists_admissibleScales` / `DiagonalScale.exists_diagonal_scales`
finite-prefix scale selection once normalized-template constants `C[j,m]` are
supplied.  Its remaining analytic ownership seam is that a plain callable can
supply those numbers without recording which hierarchy coefficients justify
them.

This increment closes only that provenance/ownership seam.  It introduces a
dependency-tracked finite-prefix certificate in
`background_hierarchy_cutoff_bounds.py`: every `C[j,m]` must name one hierarchy
identity, one source revision, and a nonempty coefficient dependency chain that
includes coefficient order `j` and contains no future order.  A certificate is
accepted only when it contains the exact triangular set

`1 <= j <= J`, `0 <= m <= j+2`.

Only after those checks may the bounds enter the existing SlowBorel cutoff
scheduler.

The implementation baseline is repository `main`
`1d812ff01d4575551f804c3599d6efa8985c7d2f`, which includes the hierarchy-owned
angular and axial third-eta preceding-diffusion bridges.  Open Agent-2 work on
stronger sixth/fifth mixed hierarchy jets is intentionally not duplicated.

## Mathematical boundary

This does **not** derive the paper's uniform `C[j,m]` constants from the current
finite hierarchy.  Instead it provides the fail-closed contract that a future
analytic hierarchy-bound layer must satisfy.  In particular:

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
provider fixture, not paper coefficient data.  It checks:

1. the provider is queried for exactly the triangular SlowBorel index set;
2. the bound rows passed to the scheduler are exactly the certified rows;
3. every scheduled dyadic edge satisfies the existing log-margin check;
4. every accepted bound carries a same-order coefficient dependency;
5. missing/duplicate indices, future dependencies, missing own-order
   dependencies, nonfinite/nonpositive values, hierarchy/revision mismatches,
   and wrong provider indices fail closed; and
6. `paper_exact` remains false at both certificate and schedule-witness levels.

A local isolated regression of this new module/test payload passed 9 tests
before publication; repository CI remains authoritative for integration with
the real scheduler and complete suite.

## Remaining boundary for Agent 2 / Agent 7

Stage 2 remains `formal-structure`, with `full_reconstruction=false` and
`paper_exact_velocity_available=false`.

The next coefficient-level requirement for Agent 2 is still the genuine
hierarchy-owned derivative data needed to finish `actualLowerSource'''`, then
`f_n'''` / higher Picard rows and actual hierarchy-derived coefficients.  In
parallel, Agent 7 still needs a real analytic provider that *proves* uniform
`C[j,m]` bounds from those hierarchy coefficients, rather than fixtures or
caller assertions.  Only after that provider exists can the current certificate
be populated by paper-derived data and extended from finite prefixes to the
infinite recursive cutoff/all-jets-flat argument, retained recurrence
cancellation identities, first-omitted-order residual factorization, and
Proposition 5.3 arbitrary-order/super-algebraic decay.

A green test for this finite-prefix bridge is therefore not evidence of
all-order convergence or a paper-exact background.
