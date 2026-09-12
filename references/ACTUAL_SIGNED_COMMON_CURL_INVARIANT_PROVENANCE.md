# Actual signed invariant-to-common-curl jets provenance

Status: **formal-structure**. This increment does not make the reconstructed velocity paper-exact.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Jets file: `NavierStokes/ActualSignedStageControls.lean`
- Jets theorem: `ActualSignedStageControls.fullRequest_jets_from_invariant`
- Curl file: `NavierStokes/ActualSignedCommonDynamics.lean`
- Curl theorem: `ActualSignedCommonDynamics.common_curl_and_divergence`

The pinned jets theorem starts from `CorrectionStep.CycleAnalyticInvariant` and proves `PeriodizedWaveBounds.UniformLocalJets` for the literal `LocalSignedRequest.fullRequest`. Its jet domain is an arbitrary `K`; for the common-curl application the relevant instantiation is exactly `ActualSignedStageControls.phaseCell`. The theorem returns exponent `sigma` and strip `HarmonicWaveInteraction.productStrip G.strip`.

The actual common-curl theorem consumes `UniformLocalJets` on `ActualSignedStageControls.fullStrip`, where the pinned source defines `fullStrip := ActualPrimaryBounds.fullStrip`. Therefore a valid composition must additionally carry theorem-grade evidence that the `SignedMeanGain.Geometry` strip is the actual primary strip, so the invariant-derived jet statement is the same full-strip statement consumed by `common_curl_and_divergence`.

## Repository increment

`src/openai_ns_reconstruction/actual_signed_common_curl_invariant.py` adds a fail-closed composition gate between the already-landed `ActualSignedCommonCurlAdmission` and an external theorem application of `fullRequest_jets_from_invariant`.

The gate requires exact agreement on `B`, `N0`, the opaque full-request identity, and the exponent identity. It also requires theorem-certified facts for the cycle analytic invariant, the literal `LocalSignedRequest.fullRequest` identity, actual-strip equality, the `phaseCell` instantiation, the invariant-sigma/exponent identity, and the jets theorem application itself. All source files, declarations, and dependency symbols are pinned exactly. Sampled, fitted, and numeric-scan evidence is rejected.

This closes one specific logical seam in the previous curl admission: its `UniformLocalJets` premise can no longer be represented as an unrelated theorem oracle when using this composition gate; it must be identified with the output of the actual cycle-invariant theorem for the same request and exponent.

## Regression boundary

`tests/test_actual_signed_common_curl_invariant.py` uses **synthetic theorem metadata only**. It verifies the valid composition and rejects cross-wired `B/N0`, request or exponent identities, missing theorem facts, stale theorem/dependency pins, and sampled/fitted/numerical evidence.

## Explicit non-claims

This increment does **not** replay Lean, authenticate an external theorem producer, verify an actual `CycleAnalyticInvariant` value, materialize the actual request or signed wave, link the signed label to a machine-exported canonical slow-box scope, or construct the final velocity. Consequently the following remain false:

- `full_request_jets_theorem_machine_replayed`
- `actual_cycle_invariant_machine_verified`
- `common_curl_theorem_machine_replayed`
- `actual_signed_wave_values_materialized`
- `paper_exact_divergence_free_wave_available`
- `paper_exact_velocity_available`

No surrogate background, profile, request, amplitude, or oscillatory wave is introduced.
