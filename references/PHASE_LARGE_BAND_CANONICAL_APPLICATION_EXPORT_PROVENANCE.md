# Section 6→7 canonical Prepared application-export handoff provenance

Status: **formal-structure**. This increment defines a fail-closed, content-addressed handoff for a future formal export of the actual `PrimaryGeometryAssembly.canonicalPrepared` application. It does not itself run Lean, enumerate the manuscript active set, materialize physical fields, or upgrade Eqs. (7.9)–(7.11).

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant symbols:

- `PrimaryGeometryAssembly.canonicalPrepared`
- `PrimaryGeometryAssembly.canonicalPhases`
- `LeadingStressWeights.FullTrueCone`

At the pinned revision, `canonicalPrepared hcone B r0 N0` returns a single `Prepared` object whose internal natural `N` is chosen by the formal construction. The repository previously checked that all admitted active boxes refer to one shared opaque Prepared identity and one `(B,r0,N0)` tuple, but it had no strict interchange contract for the future external formal artifact that should supply the actual selected object.

## Landed capability

`src/openai_ns_reconstruction/phase_large_band_canonical_application_export.py` adds a deterministic handoff contract. A supplied export must:

- declare schema `openai-ns/canonical-prepared-application-export/v1` and producer kind `lean-formal-export`;
- pin the exact formal repository and commit;
- reproduce the existing `(source_id, source_revision, prepared_instance_id)` identity;
- reproduce the canonical `(B,r0,N0)` tuple, with `r0` serialized by canonical `float.hex()` rather than decimal text;
- reproduce exactly the current canonical scope's unsigned slow-box set and signed ± family set;
- carry exactly the required canonical Prepared / phases / `FullTrueCone` symbols;
- carry a SHA-256 of the exact repository-side canonical scope contract; and
- carry a second SHA-256 over its own canonical JSON-safe payload, so any post-export metadata mutation is rejected.

The external record also carries `prepared_N`. Python only checks that it is a natural number and that it is included in the content digest. The current repository has no independent way to prove that value is the `N` selected by Lean, so no truth flag is promoted from its presence.

This design intentionally does **not** expose a production helper that fabricates a `lean-formal-export` record from the existing metadata. A future formal runner must supply the record. SHA-256 here is an integrity mechanism, not an authenticity or theorem-proof mechanism.

## Regression coverage

`tests/test_phase_large_band_canonical_application_export.py` uses the same kind of explicitly synthetic theorem metadata used by the earlier formal-structure gates. It checks:

- successful exact-scope content-addressed admission;
- rejection after changing `prepared_N` without recomputing the payload digest;
- rejection of a scope digest from another scope;
- rejection of missing unsigned boxes or signed labels even when the record is internally rehashed;
- rejection of Prepared identity or canonical-parameter drift;
- rejection of missing pinned application symbols;
- rejection of numeric-scan producer records; and
- preservation of all paper-exact / theorem-replay truth flags as false, including when an opaque `prepared_N` is present.

The test export is synthetic and explicitly says that it is not a replayed Lean artifact.

## Explicit non-claims / next blocker

A successful `LargeBandCanonicalApplicationExportAdmission` proves only that an externally supplied record is internally untampered and exactly matches the repository's admitted canonical scope. It does **not** establish who produced the record or replay its theorem terms.

The following therefore remain false:

- `actual_canonical_prepared_application_machine_verified`
- `actual_cell_index_enumeration_machine_verified`
- `actual_physical_base_values_materialized`
- `actual_base_fields_verified`
- `theorem_applications_machine_replayed`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

The next truth-raising increment requires the formal runner itself: machine-linked Lean export/replay of `canonicalPrepared`, its chosen `N`, the true active `Index W N`/`CellIndex` slice, and then field/theorem witnesses for the phase/base/frame/damping chain. No surrogate profile or hand-tuned oscillatory wave is introduced here.
