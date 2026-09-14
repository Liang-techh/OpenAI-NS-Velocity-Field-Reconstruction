# Sections 6–7 Prepared LocalBase binding provenance

## Baseline

- Repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction`
- Base `main`: `51e7b91bf5bdd9dfd52ecc74f3e1b8cb61def464`
- Scope: Issue #3, Sections 6–7 only. Sections 8–9 mean corrections/residual iteration remain Agent 8 scope and are not modified here.
- Current-main reconciliation: Agent 8's merged Section 9 finite-prefix/all-order endpoint work remains untouched. Current main additionally contains PR #300's Stage 1 actual-schedule angular reference eta jets and PR #301's Stage 2 hierarchy-owned third-eta Omega/X work; neither overlaps this four-file Sections 6–7 increment.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/PrimaryGeometryAssembly.lean`
- Field: `PrimaryGeometryAssembly.Prepared.base`
- Codomain: `PhaseEstimates.LocalBaseBounds`
- Upstream LocalBase theorem: `BaseChartJets.Estimates.localBaseBounds`

The pinned `Prepared.base` field returns LocalBase bounds for the same
`PrimaryGeometryAssembly.frequency` / `axial` fields that the existing
physical-base bridge identifies with components of `FinalSlowBase.velocity`.
This increment connects that exact theorem identity to the existing
`TheoremBackedSlowBoxBaseFieldBinding`.

## Landed capability

`phase_prepared_local_base_binding.py` adds a fail-closed theorem-shaped
application witness and a composition that:

1. requires exact slow-box, source/revision, and selected `Prepared` identity;
2. reuses the `LargeBandLocalBaseAdmission` already nested in the physical-base
   binding instead of accepting an independent `M`, scale, or derivative-ratio
   payload;
3. pins the formal symbols above and rejects sampled/fitted/numeric-scan
   evidence;
4. exposes the already-certified large-band phase-error and rounded-normal
   scalar envelopes for the identical box; and
5. preserves the shared sign-free LocalBase source for the `sigma=±1` pair.

## Truth boundary

This is `formal-structure`, not a machine replay of the noncomputable Lean
application. In particular:

- `actual_prepared_base_application_machine_replayed = false`
- `actual_physical_base_values_materialized = false`
- `uniform_eq_7_9_to_7_11_verified = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

A successful bridge does **not** manufacture a numerical background, a
tangential derivative jet, a phase/frame/damping proof, an amplitude, or a wave.
The next genuine promotion requires a machine-linked export/replay of the active
`Prepared.base` applications (or equivalent actual derivative data) and the
remaining vector/frame/damping hypotheses used by the Section 7 estimates.

## Regression coverage

`tests/test_phase_prepared_local_base_binding.py` covers coherent composition,
cross-wired box/source/Prepared identities, non-theorem evidence, missing
certificate flags, formal-symbol tampering, malformed asymptotic box keys, and
wrong composition types.

The identical source/test payload previously completed GitHub Actions run #2108
successfully on an older main. This merge-forward is based on current main
`51e7b91bf5bdd9dfd52ecc74f3e1b8cb61def464`; fresh Actions on the exact resulting
head are required before merge. Prior green CI is evidence for unchanged
source/test payload only, not a substitute for current-head CI.
