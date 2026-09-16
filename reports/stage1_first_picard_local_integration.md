# Stage 1 first-Picard local integration report

Status: **locally reviewed formal-structure increment**.

This report records the local integration evidence for the Stage 1 first-Picard
increment. The root reviewed the production formulas and independently expanded
the polynomial oracle constants with exact rational arithmetic. The current
source tree now also contains the coefficientwise x1 remainder and x2 adapters;
their focused acceptance result is recorded below.

## Repository identity

- Baseline main commit: `67ab5ecd32f80225df8bce2760b943c20a7d3d7b`
- Local integration branch: `codex/stage1-complete-slow2`
- Four-file open PR 279 prerequisite head:
  `11811fb66101bf31b5be892062f6291f3091a410`
- PR 279 prerequisite files:
  - `src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_slow2_param.py`
  - `tests/test_axis_coefficient_wide_first_picard_slow2_param.py`
  - `references/AXIS_COEFFICIENT_WIDE_FIRST_PICARD_SLOW2_PARAM_PROVENANCE.md`
  - `references/provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_param.json`

## Measured evidence

The original full suite was run as:

```text
python -m pytest -q -W error
```

It exited `0` with `1684 passed` in `1404.39 s`.

The diagnostic demo exited `0`; its report recorded
`all_diagnostic_checks_passed: true` and retained
`full_reconstruction: false`.

The strict completion audit

```text
ns-reconstruct audit --require-paper-exact
```

returned native exit code `2`, as required while the construction is
incomplete.  It continues to report
`full_reconstruction: false` and
`paper_exact_velocity_available: false`.

The base focused provenance checks passed `5` tests in `0.08 s`.

The parameter plus initial aggregate focused checks passed `14` tests in
`385.63 s`. After adding direct production-path nonzero-pressure, cancellation,
and amplitude-mismatch coverage, the final aggregate file passed `9` tests in
`2.31 s`. The final `lin2` file passed `8` tests in `22.71 s`. These focused
results were returned by the assigned Luna verification workers. The original
full run predates collection of these added tests; it is not represented as a
full run of the final expanded suite.

After reconciling the companion provenance boundaries, the root ran both final
new test files together with `tests/test_provenance_manifest.py`: `22 passed in
18.62 s`, exit `0`. The strict audit was rechecked afterward and returned native
exit `2` with both reconstruction flags false.

The aggregate production test uses independently expanded rational branch
values and tests all nine channels, including nonzero pressure and cancellation.
The `lin2` production action is separately checked on explicit diagnostic
polynomials against the exact values `2`, `-2`, `4`, and `121616/1875`.
Synthetic tests exercise local algebra only; they do not certify actual-schedule
all-order convergence or global weighted coefficient-space membership.

The newer batch covers the raw angular `lin1(x1)` branch, the angular
`quad1`/`slow1` nonlinear channels, the five-channel x1 pressure source and
pressure chain, and the full coefficientwise `naturalRemainder(x1)` plus
`second_picard_jet_pair` x2 adapter. The focused remainder/x2 acceptance check
was run exactly as:

```text
python -m pytest -q tests/test_axis_coefficient_wide_first_picard_remainder.py
```

It exited `0` with `4 passed` in `2.11 s`. Source review corrected the adapter
so only the angular path applies `naturalResolvent`; the axial path uses its
source directly. No full-suite rerun followed this increment.

The generic sparse mixed-scale algebra and formal triangular solver acceptance
check was run exactly as:

```text
python -m pytest -q tests\test_axis_coefficient_formal_solver.py -W error
```

It exited `0` with `3 passed` in `1.02 s`. The check covers independent
eta-binomial polynomial and `a^4` channels, row-zero/reference truth flags, and
the n 1/n 2 sparse channels against the landed x2 representation. No full-suite
rerun followed this increment.

## Integrated scope

The increment covers the complete mixed-scale first-Picard `slow2(x1)`
aggregate and the raw axial `lin2(x1)` branch, together with the four-file PR
279 parameter prerequisite.  The aggregate keeps ordinary powers through
`Lambda^-4`, normalized pressure-linear powers through `a^2 Lambda^-3`, and
the pressure-square `a^4 Lambda^-2` term split.  The raw `lin2(x1)` branch
keeps its ordinary reference, `Lambda^-1`, and `Lambda^-2` components plus
the normalized full-pressure `a^2 Lambda^-1` component.  It does not apply the
outer `inverseL` multiplication.

The newer coefficientwise remainder adapter binds the raw `lin1`, angular
`quad1`/`slow1`, raw `lin2`, complete `slow2`, and x1 pressure branches to one
actual `x1`; it applies the outer `inverseL` and the finite natural-resolvent
recurrence on the angular side.  Its `jet_pair(n,m,eta)` keeps ordinary
`Lambda^0` through `Lambda^-5`, normalized pressure-linear `a^2 Lambda^0`
through `a^2 Lambda^-4`, and normalized pressure-square `a^4 Lambda^-3`
channels.  `second_picard_jet_pair(n,m,eta)` then forms finite coefficientwise
`x2 = referencePair + naturalRemainder(x1)/(2 Lambda)`, adding the reference
channel and shifting the remainder channels by one inverse-Lambda power.  These
are formal-structure adapters. The generic sparse mixed-scale algebra is
implemented, and the formal triangular solver now computes requested formal
coefficients directly from lower radial rows without a fixed iteration cutoff;
neither it nor the
coefficientwise adapters establish weighted-space contraction or fixed-point
convergence.

## Remaining graph boundary

The following Stage 1 graph terms remain unresolved:

- weighted-space membership and approximation-error propagation for the generic
  mixed-scale coefficient solver, together with its connection to the existing
  `actual_schedule_picard_certificate` scalar contraction gate;
- weighted-space norm propagation and convergence certification for those
  coefficient iterates (the scalar one-step and geometric-tail bounds already
  exist in `axis_fixed_point_picard.py`);
- global weighted `AxisSpace` certification and final `phi/u` fields.

The generic mixed-scale algebra/triangular formal solver is now implemented and
its focused acceptance result is recorded above. It remains a formal local
recursion and does not provide weighted-space membership or a converged fixed
point.

The result does not claim a complete reconstruction or fixed point.  No
interval or Lean verification, remote CI result, or publication was performed
for this local integration.

All numerical checks in this report are local executable evidence and do not
promote the repository's declared `paper_exact_velocity_available` or
`full_reconstruction` flags.
