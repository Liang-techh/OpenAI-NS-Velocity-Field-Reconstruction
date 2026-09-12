# Actual signed finite-head mean-defect admission provenance

## Scope

This increment advances Issue #3 Section 8 without manufacturing a compact mean
repair. It binds the already-admitted canonical signed Section 6/7 scope to the
pinned formal theorem
`NavierStokes.ActualSignedMeanBinding.requested_cross_defect`.

At formal source commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, that theorem states on the actual
strip, for every band `n` and component `i : Fin 2`, that the actual signed mean
cross defect is exactly

`-missingWeight (choice B N0).prepared.N (physicalScale n x) * requestedStress`.

Its proof uses `ActualSignedMeanBinding.requested_cross_factor` together with
`ActualPrimaryCovariance.partitionFactor_eq_one_sub_missing`. The theorem itself
is all-band; this adapter intentionally admits only the finite head
`n <= Prepared.N`, complementary to the separately admitted exact-cancellation
tail `Prepared.N + 1 <= n`.

## Admission boundary

`actual_signed_mean_defect.py` accepts only theorem-grade metadata and requires:

- the same canonical `B/N0` and selected `Prepared.N` already admitted by the
  canonical signed scope;
- an explicitly certified finite-head band `n <= Prepared.N`;
- exact pins for `ActualSignedMeanBinding.lean`, `ActualPrimaryCovariance.lean`,
  `requested_cross_defect`, `requested_cross_factor`, and
  `partitionFactor_eq_one_sub_missing`;
- certified actual active-label identity, actual-strip membership, requested
  stress identity, missing-weight factor identity, and theorem application.

Lean cycle states and points remain opaque. Sampled, fitted, and numeric-scan
evidence is rejected. The tests intentionally use synthetic theorem metadata
only and do not replay Lean.

## Truth boundary

This remains `formal-structure` only. In particular:

- `requested_cross_defect_theorem_machine_replayed=false`;
- `actual_mean_cross_values_materialized=false`;
- `missing_weight_values_materialized=false`;
- `finite_head_mean_debt_materialized=false`;
- `compact_mean_correction_available=false`;
- `paper_exact_velocity_available=false`.

The exact identity creates a typed destination for a future formal export of the
real finite-head debt, but this increment does not evaluate `missingWeight`, does
not materialize `requestedStress`, and does not feed synthetic values into
`mean_rank_update.py`.
