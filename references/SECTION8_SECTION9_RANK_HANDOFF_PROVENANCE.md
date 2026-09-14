# Section 8 → Section 9 rank-input handoff provenance

Status: **formal-structure only**.

This increment closes one identity seam between the landed actual Section 8
post-temporal rank debt/local rank solve and the landed Section 9 finite-stage
estimate admission. It does not materialize a rank increment, a rank-stage
field, a physical mean family, or an iteration value.

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Files:

- `NavierStokes/ActualStageEstimates.lean`
- `NavierStokes/ActualPhysicalStageBounds.lean`
- `NavierStokes/ActualIntermediateDebtBounds.lean`
- `NavierStokes/LocalRankDefect.lean`

The exact chain used here is:

1. `ActualIntermediateDebtBounds.afterTemporal_debt_from_stepData` produces the
   full three-component `CorrectionState.debt` class for the literal
   post-temporal state.
2. `LocalRankDefect.RankGeometry.solved_rows` uses that same debt as the linear
   right-hand side of the local rank correction, while
   `RankGeometry.debt_eq_remainders` identifies the post-rank debt with the
   nonlinear remainder family.
3. `ActualStageEstimates.RunData.rank_class` derives the native rank-class
   estimate from the same `afterTemporal_debt_from_stepData` theorem.
4. `ActualStageEstimates.rankInput` passes `R.rank_class j` to
   `ActualPhysicalStageBounds.actualCycleRankInput`, which is the rank mean
   input consumed by the finite-stage estimate construction.

The new Python admission therefore requires one formal export to identify the
same correction-cycle state, post-temporal state and debt family on the Section
8 side, and the same `RunData` / `StageEstimates` identities on the Section 9
side. Cross-wiring any of these identities fails closed.

## Deliberate non-claims

This bridge does **not** replay Lean, reconstruct the Mathlib compact bump,
materialize the rank increment or rank-stage state, evaluate the nonlinear
remainder family, materialize `actualCycleRankInput`, or execute the Section 9
iteration. In particular it does not infer paper fields from opaque theorem
metadata.

The following remain false:

- `theorem_chain_machine_replayed`
- `actual_rank_increment_materialized`
- `actual_rank_stage_state_materialized`
- `actual_section9_rank_input_materialized`
- `section9_iteration_machine_materialized`
- `paper_exact_velocity_available`
- `full_reconstruction`

The next blocker is a paper-faithful materialization or machine replay of the
actual rank increment/rank-stage state and its remainder family, followed by a
materialized physical rank family suitable for the existing Section 9
finite-stage and selected-schedule chain. The Agent 3 wave lane remains an
upstream dependency for the actual oscillatory/correction-cycle data; no
surrogate wave or synthetic defect is promoted here.
