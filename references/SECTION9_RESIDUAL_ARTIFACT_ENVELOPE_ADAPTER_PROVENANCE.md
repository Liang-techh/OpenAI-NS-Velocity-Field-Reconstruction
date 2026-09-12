# Section 9 residual-artifact envelope adapter provenance

Status: **formal-structure only**. This file records the narrow Issue #4 increment in `section9_residual_artifact_envelope.py`; it is not evidence that the genuine Section 9 correction sequence or the paper's Eq. (9.21) residual has been reconstructed.

## Paper location and intended bridge

The adapter sits between the Section 9 residual estimate in Lemma 9.8 / Eq. (9.18), the locally finite field construction in Eq. (9.21), and the already-landed Section 9 -> Section 10 endpoint-majorant path used for the smooth extension through `t=1`. The repository source pin for the manuscript remains `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`; no new PDF hash or Lean-build claim is made here.

Integration baseline for this increment: `main` commit `3fff875b1ce9d01cf355fdd72f7a2a1bfd0b4edd` (PR #170 already integrated). The adapter reuses the existing exact-byte SHA-256 identity function in `section9_residual_artifact_fingerprint.py` and the analytic uniform-envelope conversion in `section9_endpoint_majorant_adapter.py`.

## What is newly executable

`derive_section9_endpoint_majorants_from_residual_artifact(residual_artifact, theorem_manifest)` requires two separate byte payloads. `residual_artifact` is the content-addressed object. `theorem_manifest` is a strict UTF-8 JSON sidecar whose declared `residual_artifact_sha256` must equal the SHA-256 of those exact residual bytes; the manifest itself is never substituted for the residual object.

The schema `section9-eq9.21-uniform-envelope-manifest-v1` carries one source id/revision, one Section 9 stage and spatial window, one exact rational `valid_q_upper`, and a contiguous endpoint-degree prefix `0..N`. Per degree it carries exact-rational strings for `h` and `K_m` (`derivative_loss`), integer `P_{j,m}` (`log_power`), certified leading/flat-remainder bounds, and the already-required theorem certification booleans. Unknown object fields, duplicate JSON keys, noncontiguous degrees, and residual-digest mismatches fail closed.

After parsing, every row is materialized as a `Section9UniformResidualEnvelopeWitness`, content-addressed using the exact residual bytes, and passed through the existing analytic Eq. (9.18) -> endpoint-majorant adapter. Thus the endpoint coefficient is derived by the existing stationary-point calculation from the artifact-bound theorem data instead of being supplied as a separate caller value.

## Independent regression boundary

`tests/test_section9_residual_artifact_envelope.py` uses an explicitly synthetic residual byte string and synthetic theorem-data manifest. It checks that the content digest is the residual digest rather than the manifest digest, that a one-byte residual mutation breaks the binding, that degree gaps fail, and that unknown/duplicate JSON keys fail. These tests validate the adapter contract only; synthetic fixtures are not paper evidence.

## Unresolved truth boundary

This adapter does **not** derive `C_{j,m}`, `K_m`, `P_{j,m}`, or flat-remainder bounds from the residual algebra. It does not prove that the supplied residual bytes are the manuscript's actual Eq. (9.21) residual, construct the infinite locally finite Section 9 sum, prove Proposition 9.9, produce locally uniform endpoint limits, establish all-order Borel smoothness, or construct the final smooth compact forcing. Consequently `source_majorants_derived_from_actual_residual_verified=false`, `actual_section9_sequence_verified=false`, and `paper_exact_velocity_available=false` remain mandatory.

The next substantive step is upstream of this adapter: materialize the genuine Eq. (9.21) residual/provider from the actual Section 9 sequence and generate the theorem manifest from that provider's analytic bound derivation, rather than hand-authoring its `C/K/P` and flat-remainder rows.
