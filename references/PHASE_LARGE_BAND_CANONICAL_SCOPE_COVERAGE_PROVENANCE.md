# Section 6→7 canonical Prepared scope-coverage provenance

Status: **formal-structure**. This increment does not enumerate the manuscript `CellIndex`, evaluate a Lean `canonicalPrepared`, replay theorem applications, or upgrade Eqs. (7.9)–(7.11).

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

File: `NavierStokes/PrimaryGeometryAssembly.lean`

The relevant pinned structure is:

- `PrimaryGeometryAssembly.canonicalPrepared`
- `PrimaryGeometryAssembly.canonicalPhases`
- `PrimaryGeometryAssembly.Prepared`
- `PrimaryGeometryAssembly.Index`

At the pinned revision, `canonicalPrepared hcone B r0 N0` returns one `Prepared H v (2 * activeRight W) B r0 N0`. That single `Prepared` record contains `base : ∀ L : Index W N, LocalBaseBounds ...`; `canonicalPhases` then constructs both `Fin 2` branches from the same selected object. Therefore a family-wide canonical claim requires one Prepared identity across the active boxes, not merely a collection of individually valid per-box canonical metadata records.

## Landed capability

`src/openai_ns_reconstruction/phase_large_band_canonical_scope_coverage.py` composes the existing exact physical-base scope coverage with the per-box canonical physical-base bindings from `phase_large_band_canonical_physical_base.py`.

For one admitted fixed-band scope it now requires:

- exactly one canonical binding for every unsigned slow box already present in `LargeBandPhysicalBaseScopeCoverage`;
- no missing, extra, or duplicate canonical boxes;
- each canonical binding to wrap the exact physical-base binding already admitted for that box;
- all boxes to share one `(source_id, source_revision, prepared_instance_id)` Prepared identity;
- all boxes to share one canonical constructor tuple `(B, r0, N0)`; and
- the union of all per-box ± sign pairs to equal the already admitted signed family scope.

This is a bookkeeping/integrity certificate over theorem-facing evidence. It does not infer the shared Prepared identity from sampled field values and does not create a background or oscillatory-wave surrogate.

## Independent regression coverage

`tests/test_phase_large_band_canonical_scope_coverage.py` uses two synthetic slow boxes so the global-object requirement is non-vacuous. It checks:

- successful exact-set coverage with one shared canonical Prepared identity;
- rejection of distinct per-box Prepared instance IDs even when every per-box canonical binding is individually valid;
- rejection of canonical `(B,r0,N0)` drift across boxes;
- rejection of missing or duplicate canonical box coverage;
- rejection of a canonical binding built from a different physical-base record for the same box; and
- preservation of all paper-exact truth flags as false.

The fixtures are explicitly synthetic theorem metadata. They are not manuscript `CellIndex` data, evaluated `FinalSlowBase.velocity`, or replayed Lean applications.

## Explicit non-claims / next blocker

The new gate verifies only coherence relative to supplied formal-theorem witnesses. It does **not** establish that the opaque shared `prepared_instance_id` is the value of the pinned Lean `canonicalPrepared`, and it does not enumerate the true active `Index W N`.

Therefore the following remain false:

- `actual_canonical_prepared_application_machine_verified`
- `actual_cell_index_enumeration_machine_verified`
- `actual_physical_base_values_materialized`
- `actual_base_fields_verified`
- `theorem_applications_machine_replayed`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

The next truth-raising input remains a machine-linked export/replay of the actual canonical Prepared application, including its selected `N`/active `Index W N` and field/theorem witnesses. No hand-tuned surrogate wave is introduced.
