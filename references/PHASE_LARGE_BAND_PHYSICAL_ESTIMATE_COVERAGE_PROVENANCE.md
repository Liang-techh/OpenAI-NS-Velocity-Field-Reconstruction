# Section 6→7 physical-base / estimate scope coherence provenance

Status: **formal-structure only**. This artifact composes two already landed fail-closed coverage layers on one exact theorem-facing Section 7 scope. It does not enumerate the manuscript's actual `CellIndex`, materialize a `Prepared` object or `FinalSlowBase.velocity`, replay Lean theorem applications, or make paper-exact velocity available.

## Pinned formal source reused by the composed layers

The upstream adapters remain pinned to:

- repository: `openai/NavierStokesAndEuler`
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- primary geometry module: `NavierStokes/PrimaryGeometryAssembly.lean`
- Section 7 theorem path: `BasePhaseGeometry.FamilyData.phase_estimates`, `BasePhaseGeometry.FamilyData.base_estimates`, `BasePhaseGeometry.FamilyData.coordinate_errors`, and `BasePhaseGeometry.FamilyData.damping_error`

The physical-base path already pins `BaseChartJets.CellIndex`, `PrimaryGeometryAssembly.cellDomain`, `PrimaryGeometryAssembly.family`, `PrimaryGeometryAssembly.Prepared`, the construction frequency/axial identities, the physical frequency/axial identities, and `FinalSlowBase.velocity`.

## Implemented claim

`phase_large_band_physical_estimate_coverage.py` composes:

1. `LargeBandPhysicalBaseScopeCoverage`, which requires one pinned physical-base identity binding for every unsigned slow box in an admitted primary-geometry fixed-band scope; and
2. `LargeBandFamilyEstimateCoverage`, which requires phase, coordinate-frame and damping theorem admissions to cover exactly one signed active-family scope.

The new gate mechanically requires:

- both coverage objects to carry the same complete `LargeBandActiveFamilyScopeWitness`, not merely the same number of labels;
- identical unsigned box sets and signed `sigma = ±1` label sets;
- the same `(source_id, source_revision)` on both paths;
- for each signed phase admission, resolution to exactly one sign-free physical `Prepared` binding for the same slow box;
- equality of the logical sign-free base-source box/source revision and `LocalBase` constant `M` between the estimate family and physical binding; and
- continued internal certification of both the family base-source binding and the physical-base theorem identity chain.

The returned error bounds are only the already admitted finite-scope theorem envelopes. No tolerance is fitted and no sampled field error is converted into a uniform claim.

## Why this increment is needed

Before this composition layer, physical-base exact-set coverage and phase/frame/damping exact-set coverage were individually fail-closed but independent. A caller could therefore accidentally combine coverage objects built from different theorem-facing scope witnesses and describe them as one Section 7 result. This increment rejects that cross-wiring and also checks the per-label logical base-source identity against the corresponding physical `Prepared` binding.

## Validation

Focused regression tests use synthetic theorem metadata only; they are not manuscript `CellIndex` data or exported Lean applications. They verify:

- successful one-box / two-sign composition of physical-base and all estimate paths;
- rejection of a different formal scope identity even when labels/source look the same;
- rejection when the estimate family and physical `Prepared` path use different `LocalBase M` values;
- typed, exact signed-label to physical-binding resolution; and
- preservation of every paper-exact truth gate as false.

## Explicit non-claims

A successful `LargeBandPhysicalEstimateScopeCoverage` establishes only mechanical coherence relative to the supplied theorem-facing witnesses. It does **not** independently prove the supplied fixed-band manifest is the paper's actual active `CellIndex` slice, evaluate any physical base value, or machine-replay the cited Lean theorem applications.

Therefore the following remain false:

- `actual_cell_index_enumeration_machine_verified`
- `actual_physical_base_values_materialized`
- `actual_base_fields_verified`
- `theorem_applications_machine_replayed`
- `global_all_band_uniformity_verified`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

## Remaining paper-exact boundary

The next genuine upgrade still requires a machine-linked concrete fixed-band `BaseChartJets.CellIndex` export/proof together with the corresponding `Prepared` applications or actual field witnesses. Once those real objects exist, this composition gate can require that the physical-base path and the phase/base/frame/damping theorem path refer to the same actual scope. Until then, no sampled background or surrogate oscillatory wave should be used to promote Eqs. (7.9)–(7.11).
