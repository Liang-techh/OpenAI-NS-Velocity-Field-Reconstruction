# Section 9 actual finite-stage estimates provenance

Status: **formal-structure only**.

This increment narrows the Section 9 gap between the literal correction-cycle
iteration and the already-landed residual-flatness closure admission.  It does
not construct a stage, choose a numerical decay schedule, fit a rate, or
materialize any paper field.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main file: `NavierStokes/ActualStageEstimates.lean`
- Main theorem: `NavierStokes.ActualStageEstimates.stageEstimates_of_representations`
- Ledger corollary: `NavierStokes.ActualStageEstimates.stageEstimates_ledger`
- Output structure: `NavierStokes.MixedCandidateAssembly.StageEstimates`
- Canonical finite families: `NavierStokes.ActualIterationLedger.finiteVelocity`,
  `NavierStokes.ActualIterationLedger.finitePressure`, and
  `NavierStokes.ActualIterationLedger.finiteForcing`

For `0 < q < 1`, the pinned theorem starts from one literal correction-cycle
`RunData`, its represented physical data and the checked analytic controls.  It
produces the Section 9 finite-stage `StageEstimates` for the canonical iteration
ledgers.  In particular the formal schedule is fixed by

`q1 q = q ^ (1 / 2 : Real)`,

with

`gain g = q1 q ^ (g : Real)`

and

`loss g = q1 q ^ (g : Real) * (1 + (g : Real))`.

The formal proof constructs the ledger from the literal represented
velocity/pressure/forcing and transports the estimates to the canonical
`finiteVelocity`, `finitePressure`, and `finiteForcing` families through exact
identification equalities.  The exported obligations therefore include the
open-unit-interval condition on `q`, input normalization, output strip and mean
control, stress and oscillatory-divergence control, velocity/pressure/forcing
increment estimates, and the forcing-increment margin.

## New admission boundary

`section9_actual_stage_estimates.py` accepts only a `lean-formal-export` witness.
Lean-only values are opaque stable identities.  The gate requires the same
formal application to bind the domain, forcing, `q`, literal `RunData`, physical
data, positive bump, canonical finite velocity/pressure/forcing families, the
resulting physical-data representation and `StageEstimates`, together with all
of the theorem hypotheses and exact ledger-identification certificates.

The pinned repository, commit, file, theorem/corollary/dependency symbols and
the exact `q1`/gain/loss definitions must match.  Sampled, fitted,
`numeric-scan`, or generic `formal-theorem` metadata is rejected.

## Deliberate non-claims

This admission does **not** replay the theorem in Lean, materialize an actual
`RunData` or `PhysicalData`, evaluate the finite-stage fields, or turn the
formal estimates into measured numerical bounds.  It also does not by itself
instantiate the later diagonal residual-flatness theorem.

Accordingly all of the following remain false:

- `stage_estimates_theorem_machine_replayed`
- `actual_run_data_materialized`
- `actual_physical_data_materialized`
- `finite_stage_fields_materialized`
- `stagewise_quantitative_bounds_materialized`
- `section9_iteration_machine_materialized`
- `paper_exact_velocity_available`

Regression fixtures contain synthetic formal-export metadata only.  No
surrogate profile, background, oscillatory wave, gain schedule, or fitted rate
is promoted to the paper construction.
