import hashlib
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_endpoint_majorant_adapter import (
    Section9UniformResidualEnvelopeWitness,
    derive_section9_endpoint_majorant_from_uniform_envelope,
)
from openai_ns_reconstruction.section9_endpoint_residual_source import (
    Section9ResidualEndpointSourceWitness,
)
from openai_ns_reconstruction.section9_residual_artifact_fingerprint import (
    Section9ContentAddressedResidualSourceWitness,
    Section9ContentAddressedUniformEnvelopeWitness,
    bind_content_addressed_section9_uniform_envelopes_to_residual_source,
    section9_residual_artifact_sha256,
)
from openai_ns_reconstruction.section9_stage_certificate import CertifiedBoundDatum


SOURCE_ID = "section9-eq921-residual"
SOURCE_REVISION = "fixture-revision-a"
STAGE = 20
WINDOW = 3
ARTIFACT_SHA256 = "ab" * 32
OTHER_ARTIFACT_SHA256 = "cd" * 32
ARTIFACT_BYTES = b"section9-eq921-residual-artifact\nrevision=fixture-a\n"
OTHER_ARTIFACT_BYTES = b"section9-eq921-residual-artifact\nrevision=fixture-b\n"


def _envelope(degree: int):
    return Section9UniformResidualEnvelopeWitness(
        endpoint_derivative_degree=degree,
        stage=STAGE,
        spatial_window=WINDOW,
        h=Fraction(1, 100),
        derivative_loss=Fraction(0),
        log_power=degree,
        leading_constant=CertifiedBoundDatum(
            upper_bound=2.0 + degree,
            kind="paper-derived",
            provenance=f"fixture Eq. (9.18) leading theorem degree {degree}",
        ),
        flat_remainder_uniform_bound=CertifiedBoundDatum(
            upper_bound=0.25 + 0.1 * degree,
            kind="paper-derived",
            provenance=f"fixture flat-remainder theorem degree {degree}",
        ),
        valid_q_upper=Fraction(1, 2),
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        residual_envelope_uniform_certified=True,
        flat_remainder_uniform_certified=True,
        official_support_covered_certified=True,
    )


def _source(envelopes):
    derived = tuple(
        derive_section9_endpoint_majorant_from_uniform_envelope(row)
        for row in envelopes
    )
    return Section9ResidualEndpointSourceWitness(
        stage=STAGE,
        spatial_window=WINDOW,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        evidence_kind="analytic-theorem",
        provenance="fixture theorem binding the Eq. (9.21) residual provider",
        majorant_evidence_kinds=tuple(row.bridge.evidence_kind for row in derived),
        majorant_evidence_provenance=tuple(
            row.bridge.evidence_provenance for row in derived
        ),
        past_full_jets=lambda *args: (),
        residual_provider_contract_certified=True,
        majorants_apply_to_source_certified=True,
        past_jets_same_source_certified=True,
    )


def _addressed_envelopes(envelopes, *, digest=ARTIFACT_SHA256):
    return tuple(
        Section9ContentAddressedUniformEnvelopeWitness(
            envelope=row,
            residual_artifact_sha256=digest,
        )
        for row in envelopes
    )


def test_content_addressed_binding_delegates_to_existing_source_gate():
    envelopes = tuple(_envelope(n) for n in range(3))
    addressed_source = Section9ContentAddressedResidualSourceWitness(
        source=_source(envelopes),
        residual_artifact_sha256=ARTIFACT_SHA256,
    )

    record = bind_content_addressed_section9_uniform_envelopes_to_residual_source(
        _addressed_envelopes(envelopes),
        addressed_source,
        max_endpoint_degree=2,
    )

    assert record.formal_content_addressed_binding_ready
    assert record.residual_artifact_sha256 == ARTIFACT_SHA256
    assert record.source_key == (SOURCE_ID, SOURCE_REVISION)
    assert record.binding.endpoint_degrees == (0, 1, 2)
    assert record.binding.formal_source_binding_ready
    assert record.source_majorants_derived_from_actual_residual_verified is False
    assert record.actual_section9_sequence_verified is False
    assert record.endpoint_limits_constructed is False
    assert record.all_order_borel_smoothness_verified is False
    assert record.smooth_compact_forcing_verified is False
    assert record.paper_exact_velocity_available is False


def test_same_source_revision_with_mixed_artifact_digest_fails_closed():
    envelopes = tuple(_envelope(n) for n in range(2))
    addressed = list(_addressed_envelopes(envelopes))
    addressed[1] = Section9ContentAddressedUniformEnvelopeWitness(
        envelope=envelopes[1],
        residual_artifact_sha256=OTHER_ARTIFACT_SHA256,
    )
    addressed_source = Section9ContentAddressedResidualSourceWitness(
        source=_source(envelopes),
        residual_artifact_sha256=ARTIFACT_SHA256,
    )

    with pytest.raises(ValueError, match="one content-addressed residual artifact fingerprint"):
        bind_content_addressed_section9_uniform_envelopes_to_residual_source(
            addressed,
            addressed_source,
            max_endpoint_degree=1,
        )


def test_source_digest_mismatch_fails_even_when_source_labels_match():
    envelopes = tuple(_envelope(n) for n in range(2))
    addressed_source = Section9ContentAddressedResidualSourceWitness(
        source=_source(envelopes),
        residual_artifact_sha256=OTHER_ARTIFACT_SHA256,
    )

    with pytest.raises(ValueError, match="residual_artifact_sha256 mismatch"):
        bind_content_addressed_section9_uniform_envelopes_to_residual_source(
            _addressed_envelopes(envelopes),
            addressed_source,
            max_endpoint_degree=1,
        )


def test_noncanonical_digest_is_rejected_before_binding():
    envelope = _envelope(0)
    with pytest.raises(ValueError, match="canonical lowercase 64-hex SHA-256"):
        Section9ContentAddressedUniformEnvelopeWitness(
            envelope=envelope,
            residual_artifact_sha256="A" * 64,
        )


def test_preferred_constructors_compute_digest_from_exact_raw_artifact_bytes():
    envelopes = tuple(_envelope(n) for n in range(3))
    expected_digest = hashlib.sha256(ARTIFACT_BYTES).hexdigest()
    assert section9_residual_artifact_sha256(ARTIFACT_BYTES) == expected_digest

    addressed = tuple(
        Section9ContentAddressedUniformEnvelopeWitness.from_artifact_bytes(
            envelope=row,
            residual_artifact=ARTIFACT_BYTES,
        )
        for row in envelopes
    )
    addressed_source = Section9ContentAddressedResidualSourceWitness.from_artifact_bytes(
        source=_source(envelopes),
        residual_artifact=ARTIFACT_BYTES,
    )
    record = bind_content_addressed_section9_uniform_envelopes_to_residual_source(
        addressed,
        addressed_source,
        max_endpoint_degree=2,
    )

    assert record.residual_artifact_sha256 == expected_digest
    assert all(row.residual_artifact_sha256 == expected_digest for row in addressed)
    assert addressed_source.residual_artifact_sha256 == expected_digest
    assert record.formal_content_addressed_binding_ready
    assert record.source_majorants_derived_from_actual_residual_verified is False
    assert record.actual_section9_sequence_verified is False
    assert record.paper_exact_velocity_available is False


def test_one_byte_artifact_change_fails_closed_without_manual_digest_input():
    envelopes = tuple(_envelope(n) for n in range(2))
    addressed = tuple(
        Section9ContentAddressedUniformEnvelopeWitness.from_artifact_bytes(
            envelope=row,
            residual_artifact=ARTIFACT_BYTES,
        )
        for row in envelopes
    )
    addressed_source = Section9ContentAddressedResidualSourceWitness.from_artifact_bytes(
        source=_source(envelopes),
        residual_artifact=OTHER_ARTIFACT_BYTES,
    )

    with pytest.raises(ValueError, match="residual_artifact_sha256 mismatch"):
        bind_content_addressed_section9_uniform_envelopes_to_residual_source(
            addressed,
            addressed_source,
            max_endpoint_degree=1,
        )


def test_artifact_fingerprint_rejects_text_and_empty_payloads():
    with pytest.raises(TypeError, match="exact bytes"):
        section9_residual_artifact_sha256("not raw bytes")

    with pytest.raises(ValueError, match="nonempty"):
        section9_residual_artifact_sha256(b"")
