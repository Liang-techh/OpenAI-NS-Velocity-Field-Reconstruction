# Section 9 Eq. (9.21) finite-prefix source artifact provenance

Status: **formal-structure only**.

## Scope

This increment advances Issue #4 on the real Section 9 provider path without
pretending that the genuine Eq. (9.21) residual has already been constructed.

`src/openai_ns_reconstruction/section9_eq921_prefix_source_artifact.py` consumes
the already-landed `Section9FinitePrefixEvaluationCertificate` and materializes
its admitted finite Eq. (9.21) correction-prefix data as deterministic,
content-addressed UTF-8 JSON bytes.

The artifact is intentionally and machine-readably labelled

`eq9.21-finite-prefix-source-not-residual`.

It is an upstream source artifact for a future analytic residual provider, not
a residual artifact and not an Eq. (9.18) bound certificate.

## Primary source

Paper: **Finite Time Blowup for Navier-Stokes** (OpenAI, September 2026),
Section 9 / Proposition 9.9, especially Eqs. (9.18), (9.20), and (9.21).
Canonical source pins remain in `references/SOURCES.md`.

The landed finite-prefix layer already represents the positive-stage terms

- `sum_j chi(a_j q) A_j`,
- `sum_j chi(a_j q) B_j`,
- `sum_j chi(a_j q) p_j`,

for a contiguous finite stage prefix at one `q`.  That layer explicitly does
not prove that caller-supplied stage values are the manuscript corrections and
does not construct the locally finite infinite sum.

## What this increment adds

The provider serializes exactly the finite-prefix certificate's:

- exact rational `q`, `q_big`, and each `a_j q` as integer numerator/denominator
  pairs;
- contiguous stage identities and scales;
- cutoff evaluator provenance and each stage-value provenance;
- cutoff weights, stage contributions, and accumulated `(A,B,p)` correction
  values; and
- the existing fail-closed truth boundary.

All binary64 quantities are serialized with Python `float.hex()` rather than
JSON decimal formatting.  JSON keys are sorted with compact separators, and the
raw UTF-8 bytes receive a SHA-256 content identity.  Re-materializing the same
finite-prefix certificate therefore gives the same bytes and digest; changing a
payload value changes the digest.

## Explicit non-residual boundary

A pointwise finite-prefix certificate contains no spacetime derivative jets.
It therefore cannot evaluate

`u_t + (u . grad)u - Delta u + grad p`

and cannot machine-derive the Eq. (9.18) constants/losses or the all-power flat
remainder.  The source artifact must not be supplied as if it were the genuine
residual artifact accepted by
`section9_residual_artifact_envelope.py`.

Accordingly the returned artifact always records:

- `residual_artifact_ready=false`,
- `actual_correction_field_values_verified=false`,
- `eq_9_21_infinite_sum_constructed=false`,
- `source_majorants_derived_from_actual_residual_verified=false`, and
- `paper_exact_velocity_available=false`.

This is deliberately stronger than accepting arbitrary caller-provided residual
bytes: the repository now has a canonical machine-generated identity for the
finite Eq. (9.21) source data that a future analytic residual provider must
actually differentiate.

## Regression boundary

`tests/test_section9_eq921_prefix_source_artifact.py` independently checks
deterministic byte/digest identity, exact rational serialization, binary64
hexadecimal serialization, digest sensitivity to a changed stage contribution,
and fail-closed rejection of a prefix that attempts to claim verified paper
correction values.

The fixture remains explicitly manufactured non-paper data.  These tests prove
serialization/provider behavior only.

## Remaining blocker

The next real construction step is not another digest gate.  It is to replace
the external stage-value payloads with the genuine Section 7/8 correction-field
providers plus analytic spacetime jets, construct the locally finite Eq. (9.21)
field on the common Section 9 domain, and evaluate its analytic Navier-Stokes
residual.  Only those genuine residual bytes may then drive machine-derived
Eq. (9.18) `C_{j,m}`, `K_m`, `P_{j,m}` and all-order flat-remainder bounds.

Endpoint limits, all-order Borel smoothness through `t=1`, smooth compact
forcing, finite-energy/blow-up closure, and
`paper_exact_velocity_available=true` all remain unestablished.
