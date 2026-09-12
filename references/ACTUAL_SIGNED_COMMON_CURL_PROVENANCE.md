# Actual signed common-curl theorem provenance

Status: **formal-structure**.  This record does not make the reconstructed velocity paper-exact.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/ActualSignedCommonDynamics.lean`
- Theorem: `ActualSignedCommonDynamics.common_curl_and_divergence`

The pinned theorem is stronger than the repository's earlier coefficient-level `curl_realization_algebra.py` adapter.  For an actual `SignedLabel B N0`, harmonic `n`, and point `x` in `fullStrip.domain`, under the theorem's `PeriodizedWaveBounds.UniformLocalJets` hypothesis, it proves both:

1. the cylindrical curl of the actual common signed-copy curl potential equals the corresponding `vectorMode` built from the actual `commonCorrected` amplitude; and
2. that vector mode has zero cylindrical divergence at the point.

The proof in the pinned source dispatches phase-cell points through `ActualSignedCommonDynamics.rawJets` plus `ClosedNativeWaveIdentities.native_realizes_curl_at` / `native_divergence_zero_at`, and the localized-zero branch through `LocalizedCurlRealization.native_identities_of_zero_germ`.  The common field is assembled from `ActualSignedOutputBounds.copies`, with parameters from `ActualSignedStageControls.parameters`.

## Repository increment

`src/openai_ns_reconstruction/actual_signed_common_curl.py` adds a fail-closed theorem-application admission contract.  It accepts only `analytic-theorem` or `formal-theorem` evidence and requires exact matches for the pinned repository, commit, file, theorem, and dependency-symbol tuple.  It also requires theorem-certified metadata for the theorem's uniform-local-jet hypothesis, full-strip domain membership, actual signed-copy identity, actual parameter identity, and the theorem application itself.

The Lean `SignedLabel` and `FullPoint` values are intentionally kept as opaque external identities.  Python does not invent encodings for them and does not infer them from sampled phase/background values.

## Regression boundary

`tests/test_actual_signed_common_curl.py` uses **synthetic theorem metadata only**.  It checks exact pinning, rejects sampled/fitted/numeric-scan evidence, rejects missing theorem facts, and preserves the truth boundary below.

## Explicit non-claims

The new admission does **not** run or replay Lean, authenticate an external theorem producer, materialize the actual signed wave, bind the actual signed label to the repository's current canonical slow-box scope, verify the Proposition 5.5 / Section 7 amplitude-ODE inputs, or construct the final velocity.  Consequently the following remain false:

- `theorem_application_machine_replayed`
- `actual_signed_wave_values_materialized`
- `canonical_scope_link_verified`
- `amplitude_ode_inputs_verified`
- `paper_exact_divergence_free_wave_available`
- `paper_exact_velocity_available`

No surrogate background, profile, amplitude, or oscillatory wave is introduced by this increment.
