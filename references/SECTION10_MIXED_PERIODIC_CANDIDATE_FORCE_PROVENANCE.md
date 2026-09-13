# Section 10 mixed periodic candidate-force provenance

## Scope

This record covers one bounded Issue #4 handoff: the pinned formal theorem
`NavierStokes.MixedPeriodicAssembly.exists_candidate_force` from
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
The implementation baseline for this increment is repository `main`
`051aa05d40bf490f4c0f5371542a82c73960b4fc`.

The theorem is paper-specific. It starts from the actual mixed potential/direct
velocity/pressure fields on the open past, requires smoothness, the localized
direct-field divergence condition, all-order vanishing joint jets of the
original residual, one-sided away extensions, and the axial speed blow-up path.
It then uses the fixed Section 10 spatial localization/periodization and the
paper time activation to construct a force through the locally-uniform residual
boundary-limit family.

The formal theorem returns a `CandidateProperties` witness for the activated
periodic velocity/pressure and the constructed force, proves the force globally
`C-infinity`, and identifies every force endpoint mixed jet at `t=1` with the
constructed `boundaryLimits` family.

## Why this is not the generic forcing diagnostic

`verify.reconstruct_forcing_numeric` remains a finite-difference diagnostic and
cannot certify a force by defining `f=R` and then checking `R-f=0`. The new
`section10_mixed_periodic_candidate_force.py` gate accepts only an exact
`lean-formal-export` theorem application and pins the actual formal dependency
chain (`periodicVelocity`, `boundaryLimits_locallyUniform`,
`CandidateFromLimits.force`, time activation, divergence and blow-up theorems).

## Truth boundary

This repository still does **not** replay the Lean theorem, materialize the
actual Section 7/8 -> Section 9 input fields, evaluate the locally-uniform
residual-limit family, or materialize the constructed force. Consequently this
increment does not promote any runtime artifact to a genuine final velocity or
forcing.

The following remain false:

- `actual_input_fields_materialized`
- `residual_limit_values_materialized`
- `candidate_force_field_materialized`
- `actual_section9_sequence_verified`
- `section9_field_smooth_extension_through_t1_constructed`
- `residual_artifact_ready`
- `forcing_artifact_ready`
- `divergence_free_closure_verified`
- `blow_up_closure_verified`
- `finite_energy_closure_verified`
- `endpoint_residual_closure_verified`
- `paper_exact_velocity_available`

The admission is therefore `formal-structure` only. Green CI does not change
that truth status.
