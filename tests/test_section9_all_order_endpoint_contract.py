from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_all_order_endpoint_contract import (
    SECTION10_ENDPOINT_EXACT,
    SECTION10_LATE_PLATEAU_START_EXACT,
    Section9AllOrderEndpointResidualWitness,
    admit_section9_all_order_endpoint_contract,
)


SOURCE_ID = "section9-actual-residual"
SOURCE_REVISION = "formal-export-r1"
RESIDUAL_IDENTITY = "R[u_section9,p_section9]"


def _witness(**overrides):
    kwargs = dict(
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        residual_identity=RESIDUAL_IDENTITY,
        evidence_kind="lean-formal-export",
        provenance="universal endpoint theorem fixture; no sampled evidence",
        endpoint=SECTION10_ENDPOINT_EXACT,
        late_plateau_start=SECTION10_LATE_PLATEAU_START_EXACT,
        derivative_order_quantifier="forall-natural-derivative-orders",
        analytic_majorant_template="forall n, ||partial_t D^n R|| <= M(n,t), integral M(n,t) dt < infinity",
        support_cylinder_id="section10-official-support-cylinder",
        closed_past_all_order_derivative_bounds_certified=True,
        locally_uniform_endpoint_limits_all_orders_certified=True,
        analytic_majorant_all_orders_certified=True,
        endpoint_normal_trace_match_all_orders_certified=True,
        support_exterior_zero_neighborhood_all_orders_certified=True,
        same_actual_residual_source_certified=True,
        closed_past_zero_before_certified=True,
        official_late_plateau_certified=True,
        infinite_locally_finite_borel_family_certified=True,
        max_order=None,
        finite=False,
    )
    kwargs.update(overrides)
    return Section9AllOrderEndpointResidualWitness(**kwargs)


def _admit(witness):
    return admit_section9_all_order_endpoint_contract(
        witness,
        expected_source_id=SOURCE_ID,
        expected_source_revision=SOURCE_REVISION,
        expected_residual_identity=RESIDUAL_IDENTITY,
    )


def test_universal_formal_contract_admits_without_upgrading_truth_flags():
    certificate = _admit(_witness())

    assert certificate.formal_all_order_contract_ready
    assert certificate.witness.endpoint == Fraction(1, 1)
    assert certificate.witness.late_plateau_start == Fraction(3, 4)
    assert certificate.witness.max_order is None
    assert certificate.witness.finite is False
    assert not certificate.actual_section9_sequence_verified
    assert not certificate.section9_all_order_endpoint_limits_verified
    assert not certificate.residual_artifact_ready
    assert not certificate.forcing_artifact_ready
    assert not certificate.endpoint_residual_closure_verified
    assert not certificate.paper_exact_velocity_available
    assert not certificate.full_reconstruction


def test_contract_rejects_any_finite_endpoint_ladder_claim():
    with pytest.raises(ValueError, match="max_order must be None"):
        _witness(max_order=12)

    with pytest.raises(ValueError, match="finite must be False"):
        _witness(finite=True)

    with pytest.raises(ValueError, match="universally quantify all natural orders"):
        _witness(derivative_order_quantifier="orders-0-through-100")


def test_contract_rejects_sampled_or_numerical_evidence():
    for kind in ("certified-numerical", "sampled", "fitted"):
        with pytest.raises(ValueError, match="sampled/fitted/numerical evidence is rejected"):
            _witness(evidence_kind=kind)


def test_contract_rejects_wrong_endpoint_or_late_plateau_and_float_aliases():
    with pytest.raises(ValueError, match="exact rational T=1"):
        _witness(endpoint=Fraction(999, 1000))
    with pytest.raises(ValueError, match="exact rational T=1"):
        _witness(endpoint=1.0)

    with pytest.raises(ValueError, match="exact rational 3/4"):
        _witness(late_plateau_start=Fraction(2, 3))
    with pytest.raises(ValueError, match="exact rational 3/4"):
        _witness(late_plateau_start=0.75)


def test_contract_requires_all_order_exterior_vanishing_and_trace_match():
    with pytest.raises(
        ValueError, match="support_exterior_zero_neighborhood_all_orders_certified"
    ):
        _witness(support_exterior_zero_neighborhood_all_orders_certified=False)

    with pytest.raises(ValueError, match="endpoint_normal_trace_match_all_orders_certified"):
        _witness(endpoint_normal_trace_match_all_orders_certified=False)


def test_contract_requires_actual_residual_identity_and_source_revision_match():
    witness = _witness()

    with pytest.raises(ValueError, match="source identity mismatch"):
        admit_section9_all_order_endpoint_contract(
            witness,
            expected_source_id=SOURCE_ID,
            expected_source_revision="different-revision",
            expected_residual_identity=RESIDUAL_IDENTITY,
        )

    with pytest.raises(ValueError, match="residual identity mismatch"):
        admit_section9_all_order_endpoint_contract(
            witness,
            expected_source_id=SOURCE_ID,
            expected_source_revision=SOURCE_REVISION,
            expected_residual_identity="manufactured-residual",
        )


def test_contract_requires_analytic_majorant_and_infinite_locally_finite_borel_hypotheses():
    with pytest.raises(ValueError, match="analytic_majorant_all_orders_certified"):
        _witness(analytic_majorant_all_orders_certified=False)

    with pytest.raises(ValueError, match="infinite_locally_finite_borel_family_certified"):
        _witness(infinite_locally_finite_borel_family_certified=False)
