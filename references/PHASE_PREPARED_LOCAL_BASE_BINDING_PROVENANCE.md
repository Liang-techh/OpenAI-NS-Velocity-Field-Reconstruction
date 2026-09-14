# Sections 6–7 Prepared LocalBase binding provenance

## Baseline

- Repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction`
- Base `main`: `c5d7734ce1b659d309ebe417bc88f60fc3af258e`
- Scope: Issue #3, Sections 6–7 only. Sections 8–9 mean corrections/residual iteration remain Agent 8 / parallel scope and are not modified here.
- Parallel-main delta: merged PR #298 executes the Agent 8 Section 9 finite-prefix field assembly. It is preserved unchanged and confirms that Sections 8–9 are an active parallel lane, so this increment does not touch any mean-correction/residual-iteration paths.

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

The immediately preceding head on base `1ea24db9524d34a191ebf07d3e22a7b5a1ecdd6f`
completed GitHub Actions run #2108 successfully. This rebased/replayed current-main
head requires its own fresh Actions result before merge; prior green CI is evidence
for unchanged source/test payload only, not a substitute for current-head CI.
