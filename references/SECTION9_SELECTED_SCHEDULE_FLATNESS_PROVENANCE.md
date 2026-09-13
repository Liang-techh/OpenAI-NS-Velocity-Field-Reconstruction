# Section 9 selected-schedule residual-flatness provenance

Status: **formal-structure only**.

This increment closes one formal gap after the landed actual finite-stage estimate admission. It pins the theorem in the official formal development that takes the repository's closed selected actual estimate record, chooses one common diagonal schedule, preserves the three cut bounds, and obtains every derivative order and every nonnegative power of residual flatness. It does not evaluate the chosen schedule, a stage field, a residual, or the final paper velocity.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main file: `NavierStokes/LocalResidualFlatness.lean`
- Closed theorem: `NavierStokes.LocalResidualFlatness.selected_schedule`
- General schedule theorem: `NavierStokes.LocalResidualFlatness.exists_schedule_all_jetRates`
- Quantitative conclusion: `NavierStokes.LocalResidualFlatness.AllResidualJetRates`
- Closed estimate record: `NavierStokes.LocalResidualFlatness.selectedEstimates`
- Actual estimate constructor: `NavierStokes.ActualCandidateAssembly.estimates`

At the pinned revision, `selectedEstimates` is defined from `ActualCandidateAssembly.estimates` at `ActualCandidateConstruction.selectedBudget` and `selectedThreshold`. The theorem `selected_schedule` then applies `exists_schedule_all_jetRates` to this closed record and returns one existential schedule `a : ℕ → ℕ` carrying all three of the following on the actual selected potential/direct/pressure stage families: a `MixedCandidateWitness.SelectedSchedule`, `MixedDiagonalSchedule.ThreeCutBounds`, and `AllResidualJetRates`.

`AllResidualJetRates` quantifies over every derivative order `m` and every real decay order `r ≥ 0` for the unlocalized mixed Navier--Stokes residual near the singular endpoint. The schedule is chosen outside those two quantifiers. This is materially stronger than admitting isolated stage estimates, but it remains a theorem-level formal statement rather than a materialized numerical or field reconstruction.

## New admission boundary

`section9_selected_schedule_flatness.py` accepts only a `lean-formal-export` witness for the exact closed `selected_schedule` application. Lean-only definitions and the existential schedule remain opaque stable identities. The gate requires exact identification of `selectedEstimates`, `selectedQbig`, the selected potential/direct/pressure stage families, the exported schedule witness, the three-cut-bounds conclusion, and the all-residual-jet-rates conclusion.

The pinned repository, commit, file, theorem, and dependency symbol tuple must match exactly. Sampled, fitted, `numeric-scan`, or generic `formal-theorem` metadata is rejected.

## Deliberate non-claims

This admission does **not** replay `selected_schedule` in Lean, reveal or evaluate its existential schedule, materialize any finite-stage field, produce numerical residual-jet constants, or construct the paper-exact velocity. It also does not replace the separate abstract `DiagonalResidual.allJetsFlat_residual_of_stages` admission; rather, it records the stronger closed actual-schedule theorem already present in the pinned formal source.

Accordingly all of the following remain false:

- `selected_schedule_theorem_machine_replayed`
- `selected_schedule_values_materialized`
- `finite_stage_fields_materialized`
- `residual_jet_rate_values_materialized`
- `section9_iteration_machine_materialized`
- `paper_exact_velocity_available`

Regression fixtures use synthetic formal-export metadata only. No surrogate profile, background, oscillatory wave, schedule, or fitted rate is promoted to the paper construction.
