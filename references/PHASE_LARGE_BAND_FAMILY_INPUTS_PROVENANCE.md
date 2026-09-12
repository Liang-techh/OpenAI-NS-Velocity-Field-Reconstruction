# Section 7 large-band family-input admission provenance

Status: **formal-structure**. This increment does not materialize the manuscript Section 6 partition, Proposition 5.5 background, an actual phase normal, or an oscillatory wave.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/BasePhaseGeometry.lean`
- Relevant theorem family: `BasePhaseGeometry.FamilyData.phase_estimates`, followed by `FamilyData.base_estimates`, `FamilyData.coordinate_errors`, and the damping estimate.

The pinned `FamilyData.phase_estimates` theorem assumes, before pointwise phase algebra is used, fixed family inequalities including `0 < r0`, `1 <= M`, `|u| <= M`, `1/(2*r0) <= M`, `4*r0*Tg <= M`, an active family index satisfying `LargeBand`, slow-point membership in the carrier, and slot-point membership. It then derives uniform normal and normal-motion bounds on that active carrier/slot. The later `coordinate_errors` theorem strengthens to `0 < u <= M` and combines those phase estimates with base estimates and the family identities.

## Why this layer exists

The repository already had:

1. underflow-safe large-band LocalBase admission (`phase_large_band_local_base.py`);
2. stable sign-independent base-source/revision binding (`phase_large_band_base_source.py`); and
3. the scalar large-band phase/frame/damping gate (`phase_large_band_frame.py`).

Those objects could still be combined with a detached signed label or with carrier/slot claims from an unrelated family revision. In particular, nothing typed the assertion that a specific active signed family index, its carrier membership, its slot membership, its frozen representative data, the reference scale `B`, and the viscosity coefficient all belong to the same base-source revision already used by LocalBase.

`src/openai_ns_reconstruction/phase_large_band_family_inputs.py` closes only that identity/admission gap.

## Interface and fail-closed checks

`LargeBandFamilyInputWitness` records one asymptotic signed label together with:

- exact `M`, `u`, `B`, viscosity, reference radius `r0`, and global slot constant `Tg`;
- a stable `(source_id, source_revision)`;
- theorem evidence kind and nonempty provenance;
- theorem-certified family-index/source identity;
- theorem-certified carrier and slot membership;
- theorem-certified representative-data identity;
- theorem-certified reference-scale and viscosity identity.

The constructor independently checks the scalar family prerequisites used on the route to `phase_estimates`/`coordinate_errors`: `M>=1`, `0<u<=M`, `0<B<=M`, `0<=viscosity<=4`, `r0>0`, `1/(2*r0)<=M`, and `4*r0*Tg<=M`. Only `analytic-theorem` and `formal-theorem` evidence are accepted; sampled, fitted, and numeric-scan evidence is rejected.

`LargeBandFamilyInputAdmission` then binds that witness to both an existing `LargeBandBaseSourceBinding` and an existing `LargeBandPhaseFrameCertificate`. It requires exact signed-label/box identity, exact source revision, exact `M`/band/`u`/`B`/viscosity agreement, and rechecks that both upstream gates remain intact. The exposed `family_phase_bound` is the already-landed family scalar envelope `BasePhaseGeometry.phaseConstant(M)/S_*`; it is not a measured error and does not claim that the upstream theorem facts have been machine-proved in this repository.

## Independent regression

`tests/test_phase_large_band_family_inputs.py` uses formal fixtures only. It checks that:

- either sign can be admitted while the sign-free source revision remains fixed;
- a successful admission has one exact source key and the same family scalar bound as the existing frame certificate;
- box/source/`M`/`u`/`B`/viscosity mismatches fail closed;
- the scalar global family inequalities fail at their boundaries;
- sampled/fitted/numeric-scan evidence, blank provenance/source keys, or any missing theorem fact are rejected;
- successful admission does not upgrade partition, base-field, uniform Eq. (7.9)-(7.11), or paper-exact-velocity truth flags.

The fixture values are not manuscript-flow samples and are not surrogate waves.

## Remaining boundary

This layer does **not** prove the theorem assertions stored in the witness. Still required are:

- a genuine Proposition 5.5/background provider pinned to a real source revision;
- the actual asymptotic Section 6 partition/carriers/slots and a theorem linking their active indices to that provider;
- machine-linked `FamilyData.phase_estimates` hypotheses on every active carrier/slot;
- the actual normal/slot-derivative inequalities and base comparison inequalities needed by `coordinate_errors` and damping;
- stress-cone amplitudes, curl-realized divergence-free waves, compact mean corrections, and Section 9 residual iteration built from those actual fields.

Accordingly `actual_partition_verified=false`, `actual_base_fields_verified=false`, `family_phase_estimates_machine_verified=false`, `uniform_eq_7_9_to_7_11_verified=false`, and `paper_exact_velocity_available=false` remain mandatory.
