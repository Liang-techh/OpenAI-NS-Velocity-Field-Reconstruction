import hashlib
import json

import pytest

from openai_ns_reconstruction.section9_residual_artifact_envelope import (
    derive_section9_endpoint_majorants_from_residual_artifact,
)
from openai_ns_reconstruction.section9_residual_artifact_fingerprint import (
    section9_residual_artifact_sha256,
)


def _manifest(residual: bytes, degrees=(0, 1)) -> bytes:
    rows = []
    for degree in degrees:
        rows.append(
            {
                "endpoint_derivative_degree": degree,
                "h": "1/4",
                "derivative_loss": "0",
                "log_power": degree,
                "leading_constant": {
                    "upper_bound": str(2 + degree),
                    "kind": "paper-derived",
                    "provenance": f"synthetic regression Eq. (9.18) C degree {degree}",
                },
                "flat_remainder_uniform_bound": {
                    "upper_bound": "0.125",
                    "kind": "paper-derived",
                    "provenance": f"synthetic regression flat remainder degree {degree}",
                },
                "certifications": {
                    "residual_envelope_uniform_certified": True,
                    "flat_remainder_uniform_certified": True,
                    "official_support_covered_certified": True,
                },
            }
        )
    payload = {
        "schema": "section9-eq9.21-uniform-envelope-manifest-v1",
        "residual_artifact_sha256": hashlib.sha256(residual).hexdigest(),
        "source_id": "synthetic-section9-regression",
        "source_revision": "fixture-v1",
        "stage": 20,
        "spatial_window": 3,
        "valid_q_upper": "1/2",
        "envelopes": rows,
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")


def test_exact_residual_bytes_drive_digest_envelopes_and_machine_majorants() -> None:
    residual = b"synthetic Eq. (9.21) residual artifact\x00fixture-v1"
    manifest = _manifest(residual)

    bundle = derive_section9_endpoint_majorants_from_residual_artifact(
        residual, manifest
    )

    assert bundle.formal_artifact_adapter_ready
    assert bundle.residual_artifact_sha256 == section9_residual_artifact_sha256(
        residual
    )
    assert bundle.residual_artifact_sha256 != hashlib.sha256(manifest).hexdigest()
    assert bundle.theorem_manifest_sha256 == hashlib.sha256(manifest).hexdigest()
    assert bundle.source_key == ("synthetic-section9-regression", "fixture-v1")
    assert bundle.max_endpoint_degree == 1
    assert tuple(row.endpoint_derivative_degree for row in bundle.envelopes) == (0, 1)
    assert all(
        row.residual_artifact_sha256 == bundle.residual_artifact_sha256
        for row in bundle.content_addressed_envelopes
    )
    assert all(row.formal_adapter_ready for row in bundle.derived_majorants)
    assert all(
        row.derived_uniform_coefficient > 0.0 for row in bundle.derived_majorants
    )

    assert bundle.source_majorants_derived_from_actual_residual_verified is False
    assert bundle.actual_section9_sequence_verified is False
    assert bundle.endpoint_limits_constructed is False
    assert bundle.all_order_borel_smoothness_verified is False
    assert bundle.smooth_compact_forcing_verified is False
    assert bundle.paper_exact_velocity_available is False


def test_manifest_digest_rejects_one_byte_residual_mutation() -> None:
    residual = b"synthetic-residual-v1"
    manifest = _manifest(residual)

    with pytest.raises(ValueError, match="does not match exact residual bytes"):
        derive_section9_endpoint_majorants_from_residual_artifact(
            residual + b"!", manifest
        )


def test_manifest_requires_contiguous_endpoint_degree_prefix() -> None:
    residual = b"synthetic-residual-v1"

    with pytest.raises(ValueError, match="contiguous endpoint-degree prefix"):
        derive_section9_endpoint_majorants_from_residual_artifact(
            residual, _manifest(residual, degrees=(0, 2))
        )


def test_manifest_rejects_unknown_fields_and_duplicate_json_keys() -> None:
    residual = b"synthetic-residual-v1"
    payload = json.loads(_manifest(residual).decode("utf-8"))
    payload["caller_override"] = True
    bad_extra = json.dumps(payload, sort_keys=True).encode("utf-8")

    with pytest.raises(ValueError, match="keys mismatch"):
        derive_section9_endpoint_majorants_from_residual_artifact(
            residual, bad_extra
        )

    digest = hashlib.sha256(residual).hexdigest()
    duplicate = (
        "{"
        '"schema":"section9-eq9.21-uniform-envelope-manifest-v1",'
        '"schema":"section9-eq9.21-uniform-envelope-manifest-v1",'
        f'"residual_artifact_sha256":"{digest}",'
        '"source_id":"synthetic","source_revision":"v1",'
        '"stage":20,"spatial_window":0,"valid_q_upper":"1/2","envelopes":[]'
        "}"
    ).encode("utf-8")

    with pytest.raises(ValueError, match="duplicate JSON key"):
        derive_section9_endpoint_majorants_from_residual_artifact(
            residual, duplicate
        )
