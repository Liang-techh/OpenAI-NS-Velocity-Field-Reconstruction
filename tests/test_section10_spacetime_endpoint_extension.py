from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section10_spacetime_endpoint_extension import (
    PINNED_EXTENSION_THEOREM,
    REQUIRED_FORMAL_SYMBOLS,
    Section10SpacetimeEndpointExtensionAdmission,
    Section10SpacetimeEndpointExtensionWitness,
)


def _witness(**overrides):
    data = dict(
        section9_field_id="actual-section9-field@selected-schedule",
        joint_jet_family_id="J:actual-section9-joint-jets",
        endpoint_limit_jet_id="L:actual-section9-endpoint-limits",
        closed_past_extension_id="extendTrace:T=1",
        global_extension_id="smoothExtension:T=1",
        application_id="SpacetimeGluing.exists_smooth_periodic_extension_of_limits:T=1",
        producer_kind="lean-formal-export",
        provenance="formal export fixture for the pinned endpoint theorem chain",
        endpoint_time=Fraction(1, 1),
        actual_section9_field_identification_certified=True,
        open_past_value_identity_certified=True,
        derivative_recurrence_certified=True,
        locally_uniform_left_limits_certified=True,
        unit_spatial_periods_on_open_past_certified=True,
        theorem_application_certified=True,
        global_contdiff_output_certified=True,
        agrees_with_open_past_output_certified=True,
        global_unit_spatial_periods_output_certified=True,
        right_tail_zero_from_t_plus_one_certified=True,
        endpoint_mixed_jets_equal_limits_certified=True,
    )
    data.update(overrides)
    return Section10SpacetimeEndpointExtensionWitness(**data)


def test_pinned_spacetime_endpoint_extension_admits_only_formal_structure():
    witness = _witness()
    admission = Section10SpacetimeEndpointExtensionAdmission(witness)

    assert witness.theorem_symbol == PINNED_EXTENSION_THEOREM
    assert witness.dependency_symbols == REQUIRED_FORMAL_SYMBOLS
    assert admission.endpoint_time == Fraction(1, 1)
    assert admission.certified_zero_from_time == Fraction(2, 1)
    assert admission.constructive_spacetime_extension_theorem_admitted
    assert admission.actual_left_endpoint_limit_chain_admitted
    assert admission.global_smooth_periodic_extension_output_admitted
    assert admission.endpoint_mixed_jet_preservation_admitted
    assert admission.status == "formal-structure"

    # The formal theorem application is admitted, but no Python-side field,
    # endpoint-limit values, residual or forcing is fabricated by this gate.
    assert not admission.endpoint_extension_theorem_machine_replayed
    assert not admission.section9_endpoint_limit_values_materialized
    assert not admission.global_extension_field_materialized
    assert not admission.actual_section9_sequence_verified
    assert not admission.section9_field_smooth_extension_through_t1_constructed
    assert not admission.residual_artifact_ready
    assert not admission.forcing_artifact_ready
    assert not admission.endpoint_residual_closure_verified
    assert not admission.paper_exact_velocity_available


def test_endpoint_must_be_exact_paper_time_one():
    with pytest.raises(TypeError, match="fractions.Fraction"):
        _witness(endpoint_time=1.0)
    with pytest.raises(ValueError, match="T=1"):
        _witness(endpoint_time=Fraction(3, 4))


def test_non_formal_evidence_is_rejected():
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind="numeric-scan")


def test_missing_theorem_hypothesis_or_output_certificate_fails_closed():
    for name in (
        "actual_section9_field_identification_certified",
        "open_past_value_identity_certified",
        "derivative_recurrence_certified",
        "locally_uniform_left_limits_certified",
        "unit_spatial_periods_on_open_past_certified",
        "theorem_application_certified",
        "global_contdiff_output_certified",
        "agrees_with_open_past_output_certified",
        "global_unit_spatial_periods_output_certified",
        "right_tail_zero_from_t_plus_one_certified",
        "endpoint_mixed_jets_equal_limits_certified",
    ):
        with pytest.raises(ValueError, match=name):
            _witness(**{name: False})


def test_source_or_dependency_drift_is_rejected():
    with pytest.raises(ValueError, match="formal_commit"):
        _witness(formal_commit="deadbeef")
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(theorem_symbol="NavierStokes.SpacetimeGluing.fake")
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(dependency_symbols=REQUIRED_FORMAL_SYMBOLS[:-1])


def test_identity_fields_are_not_optional_placeholders():
    with pytest.raises(ValueError, match="global_extension_id"):
        _witness(global_extension_id="   ")

    admitted = Section10SpacetimeEndpointExtensionAdmission(_witness())
    with pytest.raises(TypeError, match="witness"):
        replace(admitted, witness=object())
