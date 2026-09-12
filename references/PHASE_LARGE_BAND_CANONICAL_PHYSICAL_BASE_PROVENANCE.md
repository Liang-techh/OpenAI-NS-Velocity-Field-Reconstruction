# Section 6→7 canonical physical-base provenance

Status: **formal-structure**. This increment does not materialize a paper-exact slow base, replay a Lean application, or upgrade Eqs. (7.9)–(7.11).

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

File: `NavierStokes/PrimaryGeometryAssembly.lean`

The adapter is pinned to the following formal definitions:

- `PrimaryGeometryAssembly.canonicalPrepared`
- `PrimaryGeometryAssembly.canonicalPhases`
- `NominalConeAssembly.activeRight`
- `PrimaryGeometryAssembly.frequency`
- `PrimaryGeometryAssembly.axial`
- `FinalSlowBase.scales`
- `LeadingStressWeights.FullTrueCone`

At the pinned revision, `canonicalPrepared` removes the free box-radius choice from the generic `prepared` path by using
`upper = 2 * NominalConeAssembly.activeRight W`. `canonicalPhases` then constructs both `Fin 2` phase branches from the same canonical `Prepared` object. The `frequency` and `axial` fields used by that `Prepared.base` are built from `FinalSlowBase.scales H v upper B`.

## Landed capability

`src/openai_ns_reconstruction/phase_large_band_canonical_physical_base.py` composes the existing theorem-facing
`Prepared -> FinalSlowBase.velocity` binding with an externally certified application of the pinned canonical constructor.
It requires exact agreement of:

- `(source_id, source_revision)`;
- the opaque `prepared_instance_id` already used by the physical-base chain;
- canonical application parameters `B`, `r0>0`, and `N0`;
- the pinned repository, commit, file, and definition names; and
- separate theorem-level facts for the `FullTrueCone` input, the canonical `Prepared`/`canonicalPhases` applications, the exact `2*activeRight` upper choice, reuse of the corresponding `FinalSlowBase.scales` in both base components, and sharing of the canonical phase source across both signs.

The validator rejects sampled, fitted, and numeric-scan evidence. It does not infer any of those facts from floating-point values or from a few evaluated carrier points.

## Independent regression coverage

`tests/test_phase_large_band_canonical_physical_base.py` checks:

- successful composition with the previously landed physical-base theorem chain;
- exact source/prepared identity preservation and sign-pair reuse;
- fail-closed rejection of source revision or Prepared-instance drift;
- exact pinning of all formal symbols;
- rejection of sampled/fitted/numeric-scan evidence, blank provenance, invalid canonical parameters, and every missing certification fact; and
- preservation of all paper-exact truth flags as false.

The regression data are explicitly synthetic theorem-metadata fixtures. They are not a manuscript `CellIndex`, not a computed `canonicalPrepared` value, and not a replayed Lean proof.

## Explicit non-claims / next blocker

This increment narrows the accepted physical-base identity from an arbitrary theorem-certified `Prepared` application to the pinned canonical primary-geometry path. It still does **not** evaluate Lean's noncomputable `canonicalPrepared`, produce the actual `FullTrueCone` witness, enumerate the true `BaseChartJets.CellIndex`, materialize `FinalSlowBase.velocity`, or replay the phase/base/frame/damping theorems.

Therefore all of the following remain false:

- `actual_canonical_prepared_application_machine_verified`
- `actual_physical_base_values_materialized`
- `actual_base_fields_verified`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

The next truth-raising input remains a machine-linked formal export/replay of the actual canonical `Prepared` application and its active `CellIndex`/field witnesses. No surrogate background or hand-tuned oscillatory wave is introduced here.
