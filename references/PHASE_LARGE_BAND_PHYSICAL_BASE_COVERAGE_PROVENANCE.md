# Section 6→7 large-band physical-base scope coverage provenance

Status: **formal-structure only**. This artifact checks exact-set composition relative to an already admitted primary-geometry band slice. It does not enumerate the paper's actual `CellIndex`, materialize the slow base, or make the paper-exact velocity available.

## Upstream pinned source

This increment reuses the existing formal pins at:

- repository: `openai/NavierStokesAndEuler`
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- file: `NavierStokes/PrimaryGeometryAssembly.lean`

The upstream adapters already pin the relevant identities:

- `BaseChartJets.CellIndex`
- `PrimaryGeometryAssembly.cellDomain`
- `PrimaryGeometryAssembly.family`
- `PrimaryGeometryAssembly.Prepared`
- `PrimaryGeometryAssembly.construction_frequency`
- `PrimaryGeometryAssembly.construction_axial`
- `PrimaryGeometryAssembly.frequency_eq_physical`
- `PrimaryGeometryAssembly.axial_eq_physical`
- `FinalSlowBase.velocity`

## Implemented claim

`phase_large_band_physical_base_coverage.py` composes two previously landed fail-closed layers:

1. `LargeBandPrimaryGeometryScopeAdmission`, which binds a theorem-provenanced fixed-band unsigned box manifest to an exact signed family scope; and
2. `LargeBandPhysicalBaseBinding`, which binds one sign-free box to the pinned `Prepared` → construction frequency/axial → `FinalSlowBase.velocity` identity chain.

For one admitted band slice, the new coverage gate mechanically requires:

- exactly one physical-base binding for every unsigned box in the admitted manifest;
- no missing, extra, or duplicate box bindings;
- exact reuse of the admitted `(source_id, source_revision)` by every box;
- every per-box physical identity binding to remain internally theorem-certified; and
- the union of each binding's `sigma = ±1` pair to equal the exact signed family scope.

This prevents a strict subset of physical-base identities from being reported as family-wide coverage. The exact-set check is order-independent and contains no sampled tolerance.

## Validation

Focused regression tests use only synthetic theorem metadata fixtures. They verify:

- exact two-box coverage and sign-pair reconstruction;
- missing and extra box rejection;
- duplicate box rejection;
- source-revision drift rejection; and
- preservation of all paper-exact truth gates as false.

The fixtures are not manuscript `CellIndex` data and are not exported Lean theorem applications.

## Explicit non-claims

A successful `LargeBandPhysicalBaseScopeCoverage` proves only that Python has checked exact-set composition relative to the supplied theorem-provenanced scope. It does **not** independently verify the upstream theorem asserting that the supplied unsigned manifest is the complete paper `CellIndex` slice, evaluate any noncomputable `Prepared` object, or evaluate `FinalSlowBase.velocity`.

Therefore the following remain false:

- `actual_cell_index_enumeration_machine_verified`
- `actual_physical_base_values_materialized`
- `actual_base_fields_verified`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

## Remaining paper-exact boundary

The next genuine upgrade requires a machine-linked concrete export/proof of the fixed-band `BaseChartJets.CellIndex` slice together with the corresponding `Prepared` applications or actual field witnesses. Those witnesses can then flow through this exact-set physical-base coverage gate and the existing phase/base/frame/damping family coverage chain. Until then, no sampled background or surrogate oscillatory wave should be used to promote Eqs. (7.9)–(7.11).
