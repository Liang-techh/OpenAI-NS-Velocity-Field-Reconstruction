# Section 7 large-band `phase_estimates` admission provenance

Status: **formal-structure**. This increment does not materialize the manuscript Section 6 partition, Proposition 5.5 background, an actual phase normal/velocity field, or an oscillatory wave.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/BasePhaseGeometry.lean`
- Theorem: `BasePhaseGeometry.FamilyData.phase_estimates`

At the pinned source, `phase_estimates` states that for one active large-band family index, a slow point in its carrier, and a slot point in its slot, both

1. the actual phase normal minus `MovingFrameODE.pack (B * slope) (B • K)`, and
2. the actual `phase.velocity`

are bounded by the same `phaseConstant M / D.scale i` envelope. The proof consumes the LocalBase/representative/rounding hypotheses and the carrier/slot identities; it is not a sampled assertion.

## Why this layer exists

The repository already had an underflow-safe chain through `phase_large_band_local_base.py`, `phase_large_band_base_source.py`, `phase_large_band_frame.py`, and `phase_large_band_family_inputs.py`. That chain admitted one signed active family index and exposed the scalar `phaseConstant(M)/S_*` envelope, but it did not type a downstream theorem witness as an application of the exact pinned `FamilyData.phase_estimates` theorem for the same label and base-source revision.

`src/openai_ns_reconstruction/phase_large_band_phase_estimates.py` closes only that identity gap. A `LargeBandPhaseEstimatesWitness` must name the exact pinned repository, commit, and theorem; reuse the exact signed label and `(source_id, source_revision)` from the already-admitted family; and carry theorem-certified facts identifying the actual normal expression, actual `phase.velocity` expression, the theorem application, and uniform carrier/slot conclusion. Only `analytic-theorem` and `formal-theorem` evidence are accepted. Sampled, fitted, and numeric-scan evidence fail closed.

A successful `LargeBandPhaseEstimatesAdmission` exposes the already-landed family theorem envelope as both `normal_error_bound` and `phase_velocity_bound`. It does not measure either field and does not construct a surrogate vector.

## Independent regression

`tests/test_phase_large_band_phase_estimates.py` checks that:

- both theorem conclusions inherit exactly one common positive family envelope;
- signed-label or source-revision drift is rejected;
- stale repository/commit/theorem references are rejected;
- sampled/fitted/numeric-scan evidence and missing theorem facts are rejected; and
- successful admission leaves all field/reconstruction truth flags false.

Fixture data are formal plumbing only and are not manuscript-flow samples or surrogate waves.

## Remaining boundary

Still required before Eqs. (7.9)-(7.11) can be claimed are the genuine Section 6 active carrier/slot family and Proposition 5.5 base provider, a machine-linked proof that the admitted theorem witness actually applies to every active family index, the `base_estimates`/`coordinate_errors` and damping/frame theorem applications on those actual fields, and then stress-cone amplitudes, curl-realized divergence-free waves, compact mean corrections, and Section 9 residual improvement built from those fields.

Accordingly `actual_phase_vectors_materialized=false`, `family_phase_estimates_machine_verified=false`, `uniform_eq_7_9_to_7_11_verified=false`, and `paper_exact_velocity_available=false` remain mandatory.
