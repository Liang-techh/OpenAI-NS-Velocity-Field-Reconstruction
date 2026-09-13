# Section 8 local rank solve provenance

Status: **formal-structure only**.

This increment closes one theorem-facing seam inside Section 8: the actual
post-temporal `CorrectionState.debt` already bound to the compact rank repair is
now also bound to the pinned local variable-gauge rank solve.  It does not
materialize a correction field or infer any debt/remainder value.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main file: `NavierStokes/LocalRankDefect.lean`
- Exact five-row theorem: `NavierStokes.LocalRankDefect.RankGeometry.fiveRows`
- Linear solve theorem: `NavierStokes.LocalRankDefect.RankGeometry.solved_rows`
- Mass preservation theorem: `NavierStokes.LocalRankDefect.RankGeometry.preserve_masses`
- Post-rank debt identity: `NavierStokes.LocalRankDefect.RankGeometry.debt_eq_remainders`
- Moving support theorem: `NavierStokes.LocalRankDefect.RankGeometry.increment_supportedGauge`
- Rank geometry contract: `NavierStokes.LocalRankDefect.RankGeometry`
- Dependencies: `NavierStokes/CorrectionState.lean`, `NavierStokes/MeanRankUpdate.lean`

The crucial theorem `solved_rows` states, on the actual open slow domain and
under the explicit `RankGeometry`, local-operator, base-smoothness and slow-base
hypotheses,

`linearRows ... (rankIncrementState ...) = -CorrectionState.debt ...`.

`preserve_masses` shows the rank stage keeps the required angular and axial
moments. `debt_eq_remainders` then identifies the new rank-stage debt with the
explicit nonlinear remainder family after the linear rows have been cancelled.
`increment_supportedGauge` keeps all three rank-increment components inside the
same reserved moving support.

## Composition with current reconstruction

`section8_local_rank_solve.py` consumes an
`ActualIntermediateRankDebtAdmission`, so it can only refer to the same:

- actual correction-cycle `B/N0`;
- `ActualPrimary.commonContext B`;
- literal `ActualIntermediateDebtBounds.postTemporal x` state;
- complete `CorrectionState.debt ... (postTemporal x)` family;
- band and slow-point evaluation used by the compact rank repair;
- pinned formal repository and commit.

It additionally pins the exact `RankData`, gauge/domain identities, rank
increment/state identities and remainder-family identity. Cross-wired metadata
fails closed.

## Why this is not yet paper-exact execution

The formal theorems are pinned and their application identities are admitted,
but Python still does not materialize the noncomputable Mathlib bump, the actual
rank increment, the rank-stage state, or the nonlinear remainder values. The
existing executable `mean_rank_update.py` remains a diagnostic/formal-structure
implementation because it substitutes an explicit C-infinity bump and performs
floating matrix inversion.

The finite-head signed covariance deficit from
`ActualSignedMeanBinding.requested_cross_defect` is also not re-labelled as the
full three-component `CorrectionState.debt`; the official sources do not supply
that equality as a direct theorem, and this increment does not invent one.

Accordingly all of the following remain false:

- `theorem_chain_machine_replayed`
- `actual_rank_increment_materialized`
- `actual_rank_stage_state_materialized`
- `actual_remainder_values_materialized`
- `section9_iteration_input_materialized`
- `paper_exact_velocity_available`
- `full_reconstruction`

The next Section 8→9 blocker is to obtain a machine-linked formal export or
paper-faithful materialization of the actual rank increment/stage and its
remainder family, then feed that exact artifact into the already-landed Section
9 stage-estimate/selected-schedule chain. Regression fixtures use synthetic
opaque identities only and never count as field materialization.
