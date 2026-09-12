# Section 6→7 canonical-to-actual signed application scope provenance

Status: **formal-structure**. This increment binds the already admitted canonical `PrimaryGeometryAssembly.canonicalPrepared` signed-label roster to the already admitted `CycleAnalyticInvariant → UniformLocalJets → ActualSignedCommonDynamics.common_curl_and_divergence` application chain. It does not replay Lean, invent an encoding for `SignedLabel`, materialize a wave, or promote paper-exact velocity.

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Pinned dependency symbols:

- `PrimaryGeometryAssembly.canonicalPhases`
- `ActualSignedOutputBounds.copies`
- `ActualSignedStageControls.phaseCell`
- `LocalSignedRequest.fullRequest`
- `ActualSignedStageControls.fullRequest_jets_from_invariant`
- `ActualSignedCommonDynamics.common_curl_and_divergence`

No theorem asserting the cross-layer canonical-label → `SignedLabel` map is invented. The new link record is explicitly a future external `lean-formal-export` fact and keeps `SignedLabel` as an opaque stable identifier.

## Landed capability

`src/openai_ns_reconstruction/actual_signed_canonical_scope.py` adds three fail-closed layers:

1. `CanonicalActualSignedLabelLinkWitness` records one externally exported identity between a canonical `(ell,a,sigma)` label and one opaque Lean `SignedLabel`, pins the formal repository/commit/dependency symbols, and requires exporter certifications for canonical-family identity, phase-cell identity, and reuse of the same Prepared application.
2. `CanonicalSignedInvariantCommonCurlBinding` requires that opaque label to be exactly the one used by an already admitted `ActualSignedInvariantCommonCurlAdmission`.
3. `LargeBandCanonicalSignedCommonCurlScopeAdmission` requires exact finite signed-label coverage of the content-addressed canonical export, rejects missing/extra/duplicate canonical labels, rejects reuse of one opaque `SignedLabel` for multiple canonical labels, requires unique external application IDs, and enforces the same scope SHA-256, selected Prepared `N`, `B/N0`, repository, and commit.

The bridge deliberately consumes the existing canonical-export integrity check and the existing invariant/common-curl admission instead of reproducing their logic.

## Regression coverage

`tests/test_actual_signed_canonical_scope.py` uses synthetic theorem/formal-export metadata only. It verifies successful exact-scope composition and fail-closed rejection of missing/duplicate labels, cross-wired scope digests, wrong selected Prepared `N`, mismatched `B/N0`, opaque `SignedLabel` reuse or mismatch, sampled/fitted/numeric producer kinds, missing export facts, and dependency-symbol drift.

The fixtures are not manuscript data and are not replayed Lean artifacts.

## Explicit non-claims / next blocker

A successful scope admission proves only repository-side coherence of the handoff metadata. The following remain false:

- `actual_signed_label_export_machine_verified`
- `canonical_prepared_application_machine_replayed`
- `actual_cycle_invariant_machine_verified`
- `common_curl_theorem_machine_replayed`
- `actual_signed_wave_values_materialized`
- `paper_exact_divergence_free_wave_available`
- `paper_exact_velocity_available`

The next truth-raising step requires a real machine-linked formal export/replay that supplies the canonical Prepared application, selected `N`, true active label slice, and the actual canonical-label → `SignedLabel` identities used by the cycle-invariant/common-curl applications. This increment introduces no surrogate profile, background, amplitude, or wave.
