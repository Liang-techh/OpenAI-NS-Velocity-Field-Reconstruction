# Section 10 spacetime endpoint extension provenance

Status: **formal-structure only**. This record does not make the paper-exact velocity available.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Primary file: `NavierStokes/SpacetimeGluing.lean`
- Primary theorem: `NavierStokes.SpacetimeGluing.exists_smooth_periodic_extension_of_limits`
- Endpoint fixed here: `T = 1` exactly (`Fraction(1,1)` in Python).

The pinned theorem is the constructive spacetime result needed for the real endpoint path, not an arbitrary numerical time window. Its hypotheses are the actual open-past value identity, a full joint derivative recurrence, locally-uniform one-sided limits for every derivative, and unit spatial periodicity on the open past. Its conclusion constructs a jointly smooth global spacetime field that agrees with the original field on the open past, is unit-periodic globally, vanishes for `t >= T+1`, and preserves every full mixed endpoint jet.

The source chain used by the gate also pins:

- `SpacetimeGluing.smoothExtension`
- `SpacetimeGluing.smoothExtension_contDiff`
- `SpacetimeGluing.smoothExtension_eqOn_past`
- `SpacetimeGluing.smoothExtension_zero_from`
- `SpacetimeGluing.smoothExtension_iteratedFDeriv`
- `SpacetimeEndpoint.boundary_jets_eq_limits`
- `SpacetimeEndpoint.contDiffOn_joint_extension`
- `SpatialBorelExtension.rightExtension`
- `SpatialBorelExtension.rightExtension_right_jets`

This is stronger and more paper-specific than the already-landed `timeSwitch` support gate: the pinned `timeSwitch` is identically one near `t=1`, so it cannot itself supply the missing endpoint extension. The theorem above is the formal construction that turns genuine one-sided all-order data into a smooth global spacetime field.

## Admission rule

`Section10SpacetimeEndpointExtensionWitness` accepts only `producer_kind="lean-formal-export"`. The export must bind one actual Section 9 field identity to the theorem's joint jet family `J`, endpoint limit family `L`, the constructed closed-past extension, and the resulting global extension. It must certify every theorem hypothesis and every output property used here. Floats are rejected for the theorem endpoint; source/commit/symbol drift fails closed.

Sampled convergence, fitted rates, numeric scans, a generic `TimeWindowCutoff`, or caller-selected endpoint jets cannot pass this gate.

## Deliberate non-claims

The Python repository does **not** replay Lean here and does not materialize the exported endpoint-limit values or the global extension field. In particular this increment does not establish the missing actual Section 7/8 correction sequence, materialize the infinite/locally-finite Eq. (9.21) field, close the symmetry-axis regularity gap, build a genuine Navier-Stokes residual artifact, construct an independently supplied smooth forcing, prove endpoint residual closure, bounded energy, or the blow-up path.

Therefore the runtime truth ceiling remains fail-closed:

- `actual_section9_sequence_verified = false`
- `section9_field_smooth_extension_through_t1_constructed = false`
- `residual_artifact_ready = false`
- `forcing_artifact_ready = false`
- `endpoint_residual_closure_verified = false`
- `paper_exact_velocity_available = false`

The purpose of this increment is to replace vague "smooth extension" provenance with an exact theorem-bound handoff that a future real Section 9 exporter can satisfy without inventing a surrogate right-time branch.
