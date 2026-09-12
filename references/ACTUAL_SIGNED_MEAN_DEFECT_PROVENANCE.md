# Actual signed mean-cross defect admission provenance

## Scope

This increment advances Issue #3 Section 8 by binding the canonical signed
Section 6/7 scope to the pinned formal theorem
`NavierStokes.ActualSignedMeanBinding.requested_cross_defect`.

At formal source commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, the theorem states the exact
structural cutoff-defect identity

`meanBar(actualCross) - requestedStress = -missingWeight(Prepared.N, physicalScale) * requestedStress`.

The theorem is not a tail-only statement: given its actual-strip hypothesis it
holds for the admitted band `n`.  The separate landed
`literal_requested_cross_tail` admission handles the later fact that the
missing weight vanishes on bands satisfying `Prepared.N + 1 <= n`.  Combining
the meanings of the two formal statements explains why only a finite low-band
cutoff debt can remain, but this module does not numerically materialize that
debt.

## Pinned formal chain

The admission pins:

- `NavierStokes/ActualSignedMeanBinding.lean`;
- `ActualSignedMeanBinding.requested_cross_defect`;
- `ActualSignedMeanBinding.requested_cross_factor`; and
- `ActualPrimaryCovariance.partitionFactor_eq_one_sub_missing`.

The external application must use the same canonical `B/N0` and selected
`Prepared.N` as the already-admitted canonical signed scope.  Lean
`CorrectionState.Context`, `CorrectionState.State`, the point, and the formal
expressions for mean cross, requested stress, physical scale, and missing
weight remain opaque identifiers.  The gate accepts only
`producer_kind=lean-formal-export` and fails closed on source/symbol drift or
missing identity/application certifications.

The regression fixtures are intentionally synthetic formal-export metadata.
They verify the contract but are not a Lean replay and do not make the opaque
objects paper-exact data.

## No surrogate mean repair

No sampled/fitted/numeric-scan deficit is accepted.  The module does not assign
a numerical value to `missingWeight`, does not synthesize a stress debt from a
surrogate background, and does not feed an invented debt into any compact mean
repair.  A later compact repair may consume only a machine-linked,
theorem-derived finite-head defect export.

## Truth boundary

This increment remains `formal-structure` only:

- `requested_cross_defect_theorem_machine_replayed=false`;
- `actual_mean_cross_values_materialized=false`;
- `missing_weight_values_materialized=false`;
- `finite_head_mean_defect_values_materialized=false`;
- `finite_head_mean_defect_solved=false`;
- `compact_mean_correction_available=false`;
- `paper_exact_velocity_available=false`.
