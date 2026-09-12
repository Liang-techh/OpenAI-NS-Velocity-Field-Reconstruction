# Section 9 residual-artifact fingerprint provenance gate

## Scope

This increment hardens the already-landed Eq. (9.18) uniform-envelope to Eq. (9.21) residual-source bridge. The existing bridge freezes a semantic `(source_id, source_revision)`, Section 9 stage, spatial window, evidence kind, and provenance. Those fields remain authoritative theorem-facing metadata, but a stable revision label alone cannot detect accidental replacement of the concrete residual artifact/provider beneath that label.

`section9_residual_artifact_fingerprint.py` therefore adds a narrow content-addressed wrapper around the existing path. Every endpoint-degree Eq. (9.18) envelope is paired with a canonical lowercase 64-hex SHA-256 digest, and the frozen residual-source witness is paired with the same digest. The gate requires exact digest identity before delegating to `bind_section9_uniform_envelope_ladder_to_residual_source` for the existing degree/stage/window/evidence/provenance checks.

The preferred artifact-backed constructors now derive that digest directly from the **exact raw serialized residual-artifact bytes** via `section9_residual_artifact_sha256`. No text decoding, newline normalization, JSON reserialization, caller-written digest, or other canonicalization is used on that path. Empty payloads and non-`bytes` payloads fail closed. The lower-level digest-bearing witness constructors remain available for compatibility, but they are not the preferred path when concrete artifact bytes exist.

## What this closes

The gate fails closed when:

- different endpoint-degree envelopes name different residual-artifact digests;
- the envelope ladder and residual-source witness have different digests even if their `(source_id, source_revision)` labels are identical;
- a directly supplied digest is missing or not in canonical lowercase SHA-256 form;
- the preferred raw-byte path receives an empty or non-`bytes` artifact payload; or
- envelope and residual-source witnesses are built from bytewise different artifacts, including a one-byte mutation, without any caller-supplied digest input.

This removes the earlier preferred-path trust gap in which a caller had to type the digest independently of the artifact. It does not duplicate the existing endpoint-majorant, ladder, or residual-source validators.

## Truth boundary

The SHA-256 value remains **audit/content-identity metadata only**. Computing it from exact artifact bytes proves that the preferred wrapper is bound to those bytes; it does **not** prove that those bytes encode the manuscript's genuine Eq. (9.21) residual, that the residual was constructed from the genuine Section 7/8 fields, or that the theorem-facing Eq. (9.18) envelopes were mathematically derived from that residual.

Accordingly the following remain false/unproved:

- `source_majorants_derived_from_actual_residual_verified`;
- `actual_section9_sequence_verified`;
- the genuine Eq. (9.21) infinite correction/residual sequence and Proposition 9.9 convergence;
- endpoint-limit construction/uniqueness;
- infinite Borel right-jet convergence and all-order smoothness;
- smooth compact forcing;
- bounded-energy and blow-up closure; and
- `paper_exact_velocity_available`.

The next substantive step remains to construct or import the genuine manuscript-derived Eq. (9.21) residual artifact and bind that **same concrete byte artifact** to machine-derived uniform Eq. (9.18) data (`C_{j,m}`, `K_m`, `P_{j,m}` and all-order flat-remainder bounds). Only those independently derived bounds may flow through the content-addressed source gate toward endpoint closure.
