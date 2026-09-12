# Actual signed mean-tail admission provenance

## Scope

This increment advances Issue #3 Section 8 without manufacturing a compact mean
repair.  It binds the already-admitted canonical signed Section 6/7 scope to the
pinned formal theorem
`NavierStokes.ActualSignedMeanBinding.literal_requested_cross_tail`.

At formal source commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, that theorem states that for the
literal correction-cycle signed field, on the actual strip and tail bands
satisfying `(choice B N0).prepared.N + 1 <= n`, the appropriate mean cross is
exactly the requested stress.  Its proof passes through the actual cycle
parameters and `requested_cross_tail`; earlier normalized bands are explicitly
handled separately by `requested_cross_defect` and are not promoted to exact
cancellation here.

## Admission boundary

`actual_signed_mean_tail.py` accepts only theorem-grade metadata and requires:

- the same canonical `B/N0` and selected `Prepared.N` already admitted by the
  canonical signed scope;
- the exact tail-band threshold `Prepared.N + 1 <= n`;
- exact pins for `ActualSignedMeanBinding.lean`, `literal_requested_cross_tail`,
  `requested_cross_tail`, `cycle_cross_eq`, and `parameters_eq_cycle`;
- certified actual active-label identity, actual-strip membership, literal
  cycle-parameter identity, and theorem application.

Lean `CycleState` and point values remain opaque.  Sampled, fitted, and numeric
scan evidence is rejected.  The tests intentionally use synthetic theorem
metadata only and do not replay Lean.

## Truth boundary

This remains `formal-structure` only.  In particular:

- `literal_requested_cross_tail_theorem_machine_replayed=false`;
- `actual_mean_cross_values_materialized=false`;
- `finite_head_mean_defect_solved=false`;
- `compact_mean_correction_available=false`;
- `paper_exact_velocity_available=false`.

The finite-head cutoff deficit proved by the formal source is still real work
for Section 8.  This increment does not feed a synthetic debt into
`mean_rank_update.py`, does not identify the executable bump with Mathlib's
noncomputable bump, and does not claim a completed mean correction.
