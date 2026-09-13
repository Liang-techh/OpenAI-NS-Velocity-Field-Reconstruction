# Section 9 selected-estimates binding provenance

Status: **formal-structure only**.

This increment closes one identity seam between the landed actual finite-stage
estimate work and the landed closed selected-schedule residual-flatness
admission. It does not construct a correction stage, evaluate a field, choose
an existential diagonal schedule, or materialize the paper velocity.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Files:
  - `NavierStokes/ActualCandidateAssembly.lean`
  - `NavierStokes/LocalResidualFlatness.lean`
  - `NavierStokes/ActualStageEstimates.lean`
- Relevant symbols:
  - `NavierStokes.LocalResidualFlatness.selectedEstimates`
  - `NavierStokes.ActualCandidateAssembly.estimates`
  - `NavierStokes.GluedStageEstimates.actualStageEstimates`
  - `NavierStokes.ActualStageEstimates.stageEstimates_of_representations`
  - `NavierStokes.ActualCandidateConstruction.selectedBudget`
  - `NavierStokes.ActualCandidateConstruction.selectedThreshold`
  - `NavierStokes.ActualCandidateConstruction.selectedThreshold_geometry`
  - `NavierStokes.ActualCandidateConstruction.selectedQbig`
  - `NavierStokes.ActualCandidateAssembly.selectedPotentialStages`
  - `NavierStokes.ActualCandidateAssembly.selectedDirectStages`
  - `NavierStokes.ActualCandidateAssembly.selectedPressureStages`

At the pinned revision, `ActualCandidateAssembly.estimates B N0 hN` has output
`MixedCandidateAssembly.StageEstimates` on exactly `potentialStages B N0 hN`,
`directStages B N0 hN`, and `pressureStages B N0 hN`. Its definition delegates
to `GluedStageEstimates.actualStageEstimates`, which is the closed physical
assembly path downstream of the actual stage-estimate construction.

`LocalResidualFlatness.selectedEstimates` then specializes
`ActualCandidateAssembly.estimates` at
`ActualCandidateConstruction.selectedBudget`,
`ActualCandidateConstruction.selectedThreshold`, and
`ActualCandidateConstruction.selectedThreshold_geometry`. The same file's
`selected_schedule` theorem consumes that closed record together with
`ActualCandidateConstruction.selectedQbig` and the selected potential/direct/
pressure stage families.

## New admission boundary

`section9_selected_estimates_binding.py` requires a `lean-formal-export` witness
for that exact closed specialization and composes it with a
`Section9SelectedScheduleFlatnessAdmission`. The composition fails closed
unless the opaque identities for all five schedule inputs are identical on both
sides:

- `selectedEstimates`,
- `selectedQbig`,
- selected potential stages,
- selected direct stages,
- selected pressure stages.

The exporter must additionally certify the selected budget/threshold choice,
the selected-threshold geometry proof, the `selectedQbig` specialization, the
three selected stage-family specializations, the fact that the candidate
estimate is assembled by the actual constructor chain, and the definitional
binding from `selectedEstimates` to that closed candidate estimate. Sampled,
fitted, `numeric-scan`, and generic `formal-theorem` metadata are rejected.

This is intentionally an identity/provenance bridge rather than a second
implementation of the Section 9 estimates or schedule theorem. It prevents a
metadata record admitted by the two existing gates from silently referring to
different closed stage families or a different estimate object.

## Current-main integration note

This copy is rebuilt directly on repository `main`
`a615d64f347688d299b3d99182106817f2b8dff8`, after PR #272 landed. PR #272
adds only Section 5 exact second-eta Picard recursion infrastructure and does
not modify the Section 9 constructors used here; its landed files are preserved
unchanged. The immediately preceding PR #267 head on base
`e42cce1c499f7a2e97a303a2fc22c8d06aceb553` passed its full CI run #1897;
this current-main rebuild preserves the same bounded Section 9 implementation
and test payload while updating only the provenance baseline/integration note.

## Deliberate non-claims

This increment does **not** replay any Lean theorem or definition reduction,
materialize `selectedEstimates`, evaluate the selected stage fields, reveal the
existential schedule, produce numerical residual rates, or reconstruct the
paper velocity.

Accordingly all of the following remain false:

- `closed_estimate_binding_machine_replayed`
- `selected_estimates_values_materialized`
- `selected_stage_fields_materialized`
- `selected_schedule_values_materialized`
- `section9_iteration_machine_materialized`
- `paper_exact_velocity_available`

Regression fixtures contain synthetic formal-export metadata only. No
surrogate background, profile, wave, field, schedule, or fitted residual rate
is promoted to the paper construction.
