# Section 7 large-band `base_estimates` -> `coordinate_errors` admission provenance

Status: **formal-structure**. This increment does not materialize the manuscript Section 6 partition, Proposition 5.5 background, actual base/shear vectors, actual moving-frame coefficient errors, or an oscillatory wave.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/BasePhaseGeometry.lean`
- Theorems:
  - `BasePhaseGeometry.FamilyData.base_estimates`
  - `BasePhaseGeometry.FamilyData.coordinate_errors`

At the pinned source, `base_estimates` bounds both the base value error and the `PhaseEstimates.shearVector` error by `16 * M^2 / D.scale i` on an active carrier. The downstream `coordinate_errors` theorem consumes the previously pinned `phase_estimates` conclusion together with `base_estimates` and bounds each of `(frame i).errorA`, `(frame i).errorB`, and `(frame i).errorC` by the common `coordinateConstant M u / D.scale i` envelope on the carrier/slot.

## What this increment adds

`src/openai_ns_reconstruction/phase_large_band_coordinate_errors.py` binds those two theorem applications to the exact already-admitted `LargeBandPhaseEstimatesAdmission`.

A `LargeBandCoordinateErrorsWitness` must:

- reuse the same signed asymptotic slow label and `(source_id, source_revision)`;
- name the exact pinned Lean repository, commit, and both theorem names;
- carry theorem-certified identities for the actual base value, actual shear-vector expression, and actual frame `errorA/B/C` expressions;
- certify the `phase_estimates` and `base_estimates` dependencies of `coordinate_errors`; and
- use only `analytic-theorem` or `formal-theorem` evidence.

Sampled, fitted, and numeric-scan evidence fails closed. The admission exposes the already-transcribed scalar envelopes `16*M^2/S_*` and `coordinateConstant(M,u)/S_*`; it does not evaluate any field or substitute a surrogate vector.

## Independent regression

`tests/test_phase_large_band_coordinate_errors.py` checks that:

- the two `base_estimates` conclusions share exactly `16*M^2/S_*`;
- `errorA/errorB/errorC` inherit exactly one already-admitted frame envelope;
- signed-label or source-revision drift is rejected;
- stale repository/commit/theorem references are rejected;
- sampled evidence and every missing theorem/application identity fail closed; and
- successful admission leaves all materialization, machine-verification, uniform Eq. (7.9)-(7.11), and paper-exact truth flags false.

Fixture data are formal plumbing only and are not manuscript-flow samples or surrogate waves.

## Remaining boundary

This bridge still does **not** supply the genuine Section 6 active carrier/slot family, Proposition 5.5 paper-exact base provider, machine-linked theorem applications on all active indices, actual coordinate/frame fields, or the later damping/stress/amplitude/curl-wave construction.

Accordingly `actual_base_values_materialized=false`, `actual_shear_vectors_materialized=false`, `actual_coordinate_errors_materialized=false`, `base_estimates_machine_verified=false`, `coordinate_errors_machine_verified=false`, `uniform_eq_7_9_to_7_11_verified=false`, and `paper_exact_velocity_available=false` remain mandatory.
