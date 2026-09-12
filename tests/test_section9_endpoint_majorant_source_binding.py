from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_endpoint_majorant_adapter import (
    Section9UniformResidualEnvelopeWitness,
    derive_section9_endpoint_majorant_from_uniform_envelope,
)
from openai_ns_reconstruction.section9_endpoint_majorant_source_binding import (
    bind_section9_uniform_envelope_ladder_to_residual_source,
)
from openai_ns_reconstruction.section9_endpoint_residual_source import (
    Section9ResidualEndpointSourceWitness,
)
from openai_ns_reconstruction.section9_stage_certificate import CertifiedBoundDatum


SOURCE_ID = "section9-eq921-residual"
SOURCE_REVISION = "fixture-revision-a"
STAGE = 20
WINDOW = 3


def _envelope(degree: int, *, revision: str = SOURCE_REVISION):
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
        source_revision=revision,
        residual_envelope_uniform_certified=True,
        flat_remainder_uniform_certified=True,
        official_support_covered_certified=True,
    )


def _source(witnesses, *, revision: str = SOURCE_REVISION, evidence_kinds=None):
    derived = tuple(
        derive_section9_endpoint_majorant_from_uniform_envelope(row) for row in witnesses
    )
    kinds = tuple(row.bridge.evidence_kind for row in derived)
    if evidence_kinds is not None:
        kinds = tuple(evidence_kinds)
    return Section9ResidualEndpointSourceWitness(
        stage=STAGE,
        spatial_window=WINDOW,
        source_id=SOURCE_ID,
        source_revision=revision,
        evidence_kind="analytic-theorem",
        provenance="fixture theorem binding the Eq. (9.21) residual provider",
        majorant_evidence_kinds=kinds,
        majorant_evidence_provenance=tuple(
            row.bridge.evidence_provenance for row in derived
        ),
        past_full_jets=lambda *args: (),
        residual_provider_contract_certified=True,
        majorants_apply_to_source_certified=True,
        past_jets_same_source_certified=True,
    )


def test_uniform_envelope_ladder_is_bound_through_existing_source_gate():
    witnesses = tuple(_envelope(n) for n in range(3))
    source = _source(witnesses)

    record = bind_section9_uniform_envelope_ladder_to_residual_source(
        witnesses,
        source,
        max_endpoint_degree=2,
    )

    assert record.formal_source_binding_ready
    assert record.source_key == (SOURCE_ID, SOURCE_REVISION)
    assert record.endpoint_degrees == (0, 1, 2)
    assert record.ladder_bridge.endpoint_ladder.certified
    assert record.source_binding.formal_source_binding_ready
    assert record.source_majorants_derived_from_actual_residual_verified is False
    assert record.actual_section9_sequence_verified is False
    assert record.endpoint_limits_constructed is False
    assert record.all_order_borel_smoothness_verified is False
    assert record.smooth_compact_forcing_verified is False
    assert record.paper_exact_velocity_available is False


def test_mixed_residual_source_revisions_fail_before_ladder_binding():
    witnesses = (_envelope(0), _envelope(1, revision="fixture-revision-b"))
    source = _source((_envelope(0), _envelope(1)))

    with pytest.raises(ValueError, match="one source revision"):
        bind_section9_uniform_envelope_ladder_to_residual_source(
            witnesses,
            source,
            max_endpoint_degree=1,
        )


def test_endpoint_degree_gap_fails_closed():
    witnesses = (_envelope(0), _envelope(2))
    source = _source(witnesses)

    with pytest.raises(ValueError, match="cover exactly degrees"):
        bind_section9_uniform_envelope_ladder_to_residual_source(
            witnesses,
            source,
            max_endpoint_degree=2,
        )


def test_existing_source_gate_rejects_majorant_evidence_kind_swap():
    witnesses = tuple(_envelope(n) for n in range(2))
    source = _source(
        witnesses,
        evidence_kinds=("paper-derived", "certified-numerical"),
    )

    with pytest.raises(ValueError, match="majorant_evidence_kind_identity"):
        bind_section9_uniform_envelope_ladder_to_residual_source(
            witnesses,
            source,
            max_endpoint_degree=1,
        )
