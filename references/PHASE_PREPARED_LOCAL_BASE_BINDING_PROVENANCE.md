# Sections 6–7 Prepared LocalBase binding provenance

## Baseline

- Repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction`
- Base `main`: `065b67e9a3e8e22d697b49ccc67163db60a79783`
- Scope: Issue #3, Sections 6–7 only. Sections 8–9 mean corrections/residual iteration remain Agent 8 scope and are not modified here.
- Current-main reconciliation: the prior replay #359 was built on `1bc98e6ad4be59893db828d11b9ff0f0217ac3d8`. Current main advanced by one disjoint merged commit, PR #319, adding only the Stage 2 hierarchy-owned third-eta lower-source source/test/provenance files. None overlaps this four-file Sections 6–7 increment. Merged PR #316 and earlier Agent 8 / Stage 1/2 work remain preserved unchanged.

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

The exact preceding payload head `3702f5814c45545a8fae95eee703e493a55e40b4` completed GitHub Actions run #2311 successfully. The later replay #359 kept the same source/test payload but became stale when PR #319 advanced main from `1bc98e6ad4be59893db828d11b9ff0f0217ac3d8` to `065b67e9a3e8e22d697b49ccc67163db60a79783` while its fresh run was still queued. This replay is rebuilt directly on that latest main. Fresh Actions on the exact resulting head are authoritative before merge; prior green CI is supporting evidence only.
