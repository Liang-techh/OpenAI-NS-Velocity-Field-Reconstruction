# Section 9 residual-artifact fingerprint provenance gate

## Scope

This increment hardens the already-landed Eq. (9.18) uniform-envelope to Eq. (9.21) residual-source bridge. The existing bridge freezes a semantic `(source_id, source_revision)`, Section 9 stage, spatial window, evidence kind, and provenance. Those fields remain authoritative theorem-facing metadata, but a stable revision label alone cannot detect accidental replacement of the concrete residual artifact/provider beneath that label.

`section9_residual_artifact_fingerprint.py` therefore adds a narrow content-addressed wrapper around the existing path. Every endpoint-degree Eq. (9.18) envelope is paired with a canonical lowercase 64-hex SHA-256 digest, and the frozen residual-source witness is paired with the same digest. The gate requires exact digest identity before delegating to `bind_section9_uniform_envelope_ladder_to_residual_source` for the existing degree/stage/window/evidence/provenance checks.

## What this closes

The new gate fails closed when:

- different endpoint-degree envelopes claim different concrete residual-artifact digests;
- the envelope ladder and residual-source witness carry different digests even if their `(source_id, source_revision)` labels are identical; or
- a digest is missing or not in canonical lowercase SHA-256 form.

This prevents an integration path from silently swapping a concrete residual artifact while retaining the same human-readable source revision. It does not duplicate the existing endpoint-majorant, ladder, or residual-source validators.

## Truth boundary

The SHA-256 digest is audit metadata only. This module does **not** compute the digest from a manuscript-derived Eq. (9.21) artifact, prove that such an artifact has been constructed, or prove that the theorem-facing Eq. (9.18) envelopes were derived from it. A dishonest caller could still attach the wrong digest to theorem input; eliminating that remaining trust requires the future actual residual-construction pipeline to emit its own digest together with machine-derived bounds.

Accordingly the following remain false/unproved:

- `source_majorants_derived_from_actual_residual_verified`;
- `actual_section9_sequence_verified`;
- endpoint-limit construction/uniqueness;
- infinite Borel right-jet convergence and all-order smoothness;
- smooth compact forcing;
- bounded-energy and blow-up closure; and
- `paper_exact_velocity_available`.

The next substantive step remains unchanged: construct or import the genuine Eq. (9.21) residual revision, have that concrete artifact emit its SHA-256 identity plus machine-derived uniform Eq. (9.18) data (`C_{j,m}`, `K_m`, `P_{j,m}` and all-order flat-remainder bounds), and feed those outputs into the source-locked endpoint path.
