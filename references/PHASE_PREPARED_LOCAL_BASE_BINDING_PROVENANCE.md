# Sections 6–7 Prepared LocalBase binding provenance

## Baseline

- Repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction`
- Base `main`: `1bc98e6ad4be59893db828d11b9ff0f0217ac3d8`
- Scope: Issue #3, Sections 6–7 only. Sections 8–9 mean corrections/residual iteration remain Agent 8 scope and are not modified here.
- Current-main reconciliation: current main now also includes merged PR #316, the Section 9 actual mixed-residual all-order joint-limits bridge. The delta from the previous base `3cc438f29524acb068e40bf3a99a8293c9222628` to current main changes only the three Section 9 files for that bridge and has no overlap with this four-file Sections 6–7 increment. Existing Agent 8 PRs #309/#313/#305/#298 and Stage 1/2 work remain preserved unchanged.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/PrimaryGeometryAssembly.lean`
- Field: `PrimaryGeometryAssembly.Prepared.base`
- Codomain: `PhaseEstimates.LocalBaseBounds`
- Upstream LocalBase theorem: `BaseChartJets.Estimates.localBaseBounds`

The pinned `Prepared.base` field returns LocalBase bounds for the same `PrimaryGeometryAssembly.frequency` / `axial` fields that the existing physical-base bridge identifies with components of `FinalSlowBase.velocity`. This increment connects that theorem identity to the existing `TheoremBackedSlowBoxBaseFieldBinding`.

## Landed capability

`phase_prepared_local_base_binding.py` adds a fail-closed theorem-shaped application witness and a composition that requires exact slow-box, source/revision, and selected `Prepared` identity; reuses the `LargeBandLocalBaseAdmission` already nested in the physical-base binding rather than accepting an independent `M`, scale, or derivative-ratio payload; pins the formal symbols above and rejects sampled/fitted/numeric-scan evidence; exposes the already-certified large-band phase-error and rounded-normal scalar envelopes for the identical box; and preserves the shared sign-free LocalBase source for the `sigma=±1` pair.

## Truth boundary

This is `formal-structure`, not a machine replay of the noncomputable Lean application. In particular:

- `actual_prepared_base_application_machine_replayed = false`
- `actual_physical_base_values_materialized = false`
- `uniform_eq_7_9_to_7_11_verified = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

A successful bridge does **not** manufacture a numerical background, a tangential derivative jet, a phase/frame/damping proof, an amplitude, or a wave. The next genuine promotion requires a machine-linked export/replay of active `Prepared.base` applications (or equivalent actual derivative data) and the remaining vector/frame/damping hypotheses used by the Section 7 estimates.

## Regression and CI boundary

`tests/test_phase_prepared_local_base_binding.py` covers coherent composition, cross-wired box/source/Prepared identities, non-theorem evidence, missing certificate flags, formal-symbol tampering, malformed asymptotic box keys, and wrong composition types.

The exact preceding head `3702f5814c45545a8fae95eee703e493a55e40b4` completed GitHub Actions run #2311 successfully. Because current main later advanced through PR #316, that green run is evidence for the unchanged source/test blobs but is not substituted for current-main integration testing. This replay is built directly on current main `1bc98e6ad4be59893db828d11b9ff0f0217ac3d8`; fresh Actions on the exact resulting head are authoritative before merge.
