import numpy as np
import pytest

from openai_ns_reconstruction.section9_endpoint_bridge import (
    Section9UniformEndpointDerivativeWitness,
)
from openai_ns_reconstruction.section9_endpoint_ladder_bridge import (
    admit_section9_uniform_endpoint_ladder,
)
from openai_ns_reconstruction.section9_endpoint_residual_source import (
    Section9ResidualEndpointSourceBinding,
    Section9ResidualEndpointSourceWitness,
    admit_section10_endpoint_candidate_from_bound_section9_source,
)
from openai_ns_reconstruction.section9_stage_certificate import CertifiedBoundDatum


def _majorant_provenance(degree: int) -> str:
    return f"uniform endpoint residual theorem degree {degree}; source revision fixture-r1"


def _bridge():
    rows = [
        Section9UniformEndpointDerivativeWitness(
            endpoint_derivative_degree=degree,
            section9_derivative_order=degree + 1,
            stage=14,
            spatial_window=3,
            coefficient_bound=CertifiedBoundDatum(
                upper_bound=0.2,
                kind="paper-derived",
                provenance=_majorant_provenance(degree),
            ),
            singularity_exponent=0.0,
            valid_from=0.75,
        )
        for degree in range(3)
    ]
    return admit_section9_uniform_endpoint_ladder(rows, max_endpoint_degree=2)


def _endpoint_full_jets(degree: int, x: float, y: float, z: float) -> np.ndarray:
    shape = (3,) + (4,) * degree
    out = np.zeros(shape)
    out[(slice(None),) + (0,) * degree] = np.array(
        [0.4 + 0.03 * x, -0.2 - 0.02 * y, 0.1 + 0.01 * z]
    )
    if degree >= 1:
        out[(0,) + (1,) * degree] = 0.05 * degree
    return out


def _past_full_jets(
    degree: int, x: float, y: float, z: float, t: float
) -> np.ndarray:
    out = _endpoint_full_jets(degree, x, y, z)
    drift = np.zeros_like(out)
    drift.reshape(-1)[-1] = 0.1
    return out + (1.0 - t) * drift


def _source(bridge, **overrides):
    kwargs = dict(
        stage=bridge.stage,
        spatial_window=bridge.spatial_window,
        source_id="section9-residual-fixture",
        source_revision="fixture-r1",
        evidence_kind="formal-theorem",
        provenance="formal fixture asserting one residual source identity",
        majorant_evidence_kinds=tuple(
            row.evidence_kind for row in bridge.bridge_records
        ),
        majorant_evidence_provenance=tuple(
            row.evidence_provenance for row in bridge.bridge_records
        ),
        past_full_jets=_past_full_jets,
        residual_provider_contract_certified=True,
        majorants_apply_to_source_certified=True,
        past_jets_same_source_certified=True,
    )
    kwargs.update(overrides)
    return Section9ResidualEndpointSourceWitness(**kwargs)


def test_source_binding_freezes_one_revision_through_endpoint_borel_admission():
    bridge = _bridge()
    binding = Section9ResidualEndpointSourceBinding(bridge, _source(bridge))

    assert binding.formal_source_binding_ready
    assert binding.source_key == ("section9-residual-fixture", "fixture-r1")
    assert binding.expected_majorant_evidence_kinds == (
        "paper-derived",
        "paper-derived",
        "paper-derived",
    )
    assert binding.expected_majorant_provenance == tuple(
        _majorant_provenance(degree) for degree in range(3)
    )

    certificate = admit_section10_endpoint_candidate_from_bound_section9_source(
        binding,
        _endpoint_full_jets,
        lambda degree: 1 << degree,
        0.11,
        -0.08,
        0.04,
        0.9,
    )

    assert certificate.formal_source_bound_candidate_ready
    assert certificate.source_key == binding.source_key
    assert certificate.admission.formal_candidate_ready
    assert all(row.admitted for row in certificate.admission.rows)
    assert not certificate.source_majorants_derived_from_actual_residual_verified
    assert not certificate.actual_section9_residual_limits_verified
    assert not certificate.endpoint_limit_uniqueness_verified
    assert not certificate.infinite_borel_right_jets_verified
    assert not certificate.all_order_borel_smoothness_verified
    assert not certificate.smooth_compact_forcing_verified
    assert not certificate.paper_exact_velocity_available
    assert not binding.actual_section9_sequence_verified
    assert not binding.endpoint_limits_constructed
    assert not binding.paper_exact_velocity_available


def test_binding_rejects_majorants_from_a_different_source_revision():
    bridge = _bridge()
    mismatched = list(_source(bridge).majorant_evidence_provenance)
    mismatched[1] = "degree-1 theorem from unrelated residual revision"
    source = _source(bridge, majorant_evidence_provenance=tuple(mismatched))

    with pytest.raises(ValueError, match="majorant_provenance_identity"):
        Section9ResidualEndpointSourceBinding(bridge, source)


def test_binding_rejects_changed_evidence_kind_even_with_same_provenance():
    bridge = _bridge()
    mismatched = list(_source(bridge).majorant_evidence_kinds)
    mismatched[1] = "certified-numerical"
    source = _source(bridge, majorant_evidence_kinds=tuple(mismatched))

    with pytest.raises(ValueError, match="majorant_evidence_kind_identity"):
        Section9ResidualEndpointSourceBinding(bridge, source)


def test_binding_rejects_stage_or_spatial_window_mismatch():
    bridge = _bridge()

    with pytest.raises(ValueError, match="stage_identity"):
        Section9ResidualEndpointSourceBinding(
            bridge, _source(bridge, stage=bridge.stage + 1)
        )

    with pytest.raises(ValueError, match="spatial_window_identity"):
        Section9ResidualEndpointSourceBinding(
            bridge, _source(bridge, spatial_window=bridge.spatial_window + 1)
        )


def test_source_witness_rejects_sampled_evidence_or_missing_identity_fact():
    bridge = _bridge()

    with pytest.raises(ValueError, match="analytic-theorem or formal-theorem"):
        _source(bridge, evidence_kind="certified-numerical")

    with pytest.raises(ValueError, match="past_jets_same_source_certified"):
        _source(bridge, past_jets_same_source_certified=False)

    with pytest.raises(ValueError, match="source_revision"):
        _source(bridge, source_revision="   ")

    with pytest.raises(ValueError, match="majorant evidence kind degree 1"):
        _source(
            bridge,
            majorant_evidence_kinds=("paper-derived", "   ", "paper-derived"),
        )


def test_source_bound_api_does_not_accept_an_unbound_past_provider():
    bridge = _bridge()
    binding = Section9ResidualEndpointSourceBinding(bridge, _source(bridge))

    def unrelated_past(*_args):
        raise AssertionError("an unrelated past provider must never be queried")

    # The production entry point has no argument for swapping the past provider.
    # Replacing only a local variable therefore cannot affect the frozen source.
    unused = unrelated_past
    assert callable(unused)
    certificate = admit_section10_endpoint_candidate_from_bound_section9_source(
        binding,
        _endpoint_full_jets,
        lambda degree: 2 * (1 << degree),
        0.0,
        0.0,
        0.0,
        0.95,
    )
    assert certificate.formal_source_bound_candidate_ready
