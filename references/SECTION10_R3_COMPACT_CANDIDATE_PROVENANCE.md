# Section 10 compact whole-space candidate provenance

## Scope

This record covers one bounded Issue #4 handoff from the already-admitted
periodic Section 10 candidate to the compact whole-space model in the pinned
formalization `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
The implementation baseline for this increment is repository `main`
`30caa67ce4b9d6e95c93dd9fd0323f30298defd1`.

The pinned theorem is
`NavierStokes.R3CompactCandidate.of_localized_fields` from
`NavierStokes/R3CompactCandidate.lean`.

It consumes the same raw potential/direct/pressure fields used by the periodic
candidate.  Its periodic input is the paper time-activated candidate already
constructed by the Section 10 periodic theorem.  The theorem then keeps the
fixed spatially localized velocity and pressure before periodization, applies a
larger compact spatial cutoff to the already-constructed smooth force, and
returns the whole-space `R3CompactCandidate.Properties` witness.

Those formal output properties include smooth velocity/pressure/force, one
fixed compact spatial support for velocity and pressure, compact spatial
support for the force, zero initial velocity, compact future-time support for
the force, divergence-free velocity, the exact forced Navier--Stokes equation,
and the inherited speed-unbounded path.

## Cross-layer identity requirement

`section10_r3_compact_candidate.py` does not admit an arbitrary periodic model.
It requires a validated `Section10MixedPeriodicCandidateForceAdmission` and
checks that the raw potential/direct/pressure identities, final activated
periodic velocity/pressure identities, and theorem-constructed force identity
match that prior admission exactly.  The local-to-periodic equivalence,
support, and compact-force identities must also be certified by an exact
`lean-formal-export` theorem application.

This blocks replacing the formal periodic candidate with a manufactured field,
a generic cutoff, a fitted profile, or a residual-defined numerical force.

## Truth boundary

This is still `formal-structure` only.  Python does not replay Lean and does not
materialize the compact R3 velocity, pressure, or force.  The theorem-level
compact-support/divergence/equation/blow-up outputs are admitted as opaque
formal facts; they are not promoted to runtime-verified field artifacts.

In particular, this increment deliberately does **not** claim uniform finite
energy.  The pinned formal source contains the separate theorem
`NavierStokesR3.CompactComparisonBounds.uniformFiniteEnergy_of_compact_slab`,
which can turn a smooth field with one fixed compact support into a uniform
finite-energy bound on a closed compact time interval.  Binding that theorem to
the actual compact candidate is a separate next step.

The following remain false:

- `theorem_machine_replayed`
- `r3_compact_fields_materialized`
- `compact_force_field_materialized`
- `residual_artifact_ready`
- `forcing_artifact_ready`
- `compact_support_runtime_verified`
- `divergence_free_closure_verified`
- `blow_up_closure_verified`
- `finite_energy_closure_verified`
- `endpoint_residual_closure_verified`
- `paper_exact_velocity_available`

Green CI does not change that truth status.
