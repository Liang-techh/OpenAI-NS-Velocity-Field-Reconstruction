# Actual intermediate rank-debt provenance

Status: **formal-structure only**.

This increment narrows the remaining Section 8 gap between the actual correction
cycle and the already-landed compact rank-repair admission.  It does not create
or numerically infer any debt vector.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main file: `NavierStokes/ActualIntermediateDebtBounds.lean`
- Main theorem: `NavierStokes.ActualIntermediateDebtBounds.afterTemporal_debt_from_stepData`
- Dependency theorem: `NavierStokes.ActualIntermediateDebtBounds.stage_debt_from_stepData`
- Debt definition: `NavierStokes.CorrectionState.debt`
- Supporting file: `NavierStokes/CorrectionAnalyticStep.lean`

The pinned theorem is explicitly the full three-component slow source used by
the actual rank repair.  From checked `CorrectionAnalyticStep.StaticData`, the
same `CycleAnalyticInvariant`, actual particular inputs, and matching
`CorrectionAnalyticStep.StepData`, it proves an `UnweightedClass` bound for

`CorrectionState.debt (ActualPrimary.commonContext B) (postTemporal x)`.

The preceding `stage_debt_from_stepData` theorem derives the signed- and
temporal-stage debt bounds from the literal intermediate states rather than
accepting an estimated post-temporal debt as a premise.

## New admission boundary

`actual_intermediate_rank_debt.py` accepts only a `lean-formal-export` witness.
It keeps the Lean objects opaque while binding the following identities:

- the same `B`, `N0`, correction-cycle state, invariant and checked step data;
- the actual `ActualPrimary.commonContext B`;
- the literal `postTemporal x` state;
- the complete `CorrectionState.debt ... (postTemporal x)` family;
- one exact band/slow-point evaluation of that family.

The composition then requires the existing
`CompactMeanCorrectionAdmission` to use that exact same context, state, band,
slow point and debt instance.  Cross-wired or merely similar metadata fails
closed.

## Deliberate non-claims

This bridge does **not** identify the finite-head signed cross defect from
`ActualSignedMeanBinding.requested_cross_defect` with `CorrectionState.debt`.
Those objects occur at different layers of the construction, and no equality is
invented here.  A dedicated formal theorem/export is still required before the
finite-head defect can be called the rank debt.

It also does not materialize the post-temporal state, evaluate debt values,
replay Lean, construct the noncomputable compact bump, or produce a corrected
paper-exact velocity.

Accordingly all of the following remain false:

- `theorem_application_machine_replayed`
- `actual_cycle_state_materialized`
- `actual_post_temporal_debt_values_materialized`
- `finite_head_signed_defect_link_verified`
- `compact_mean_correction_field_materialized`
- `paper_exact_velocity_available`

Regression fixtures contain synthetic theorem metadata only.  Sampled, fitted,
and numeric-scan evidence are rejected.
