import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_natural_operator import (
    actual_schedule_chi_coefficient_state,
    actual_schedule_natural_operator,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.natural_axis import D, H
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


ROOT = Path(__file__).resolve().parents[1]


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def reference():
    return actual_schedule_reference_axis_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def natural_operator(reference):
    return actual_schedule_natural_operator(reference)


def _chi_eta(reference, eta: float) -> float:
    data = reference.reference.data
    j = reference.reference.j
    sigma = reference.reference.sigma
    h_value = H(data.h, j, eta)
    h_eta = D(data.h) + 4.0 - 12.0 * eta * eta - 2.0 * j * eta
    denominator = h_value * h_value + sigma * sigma
    return 2.0 * h_value * h_eta * sigma * sigma / (denominator * denominator)


def test_chi_state_matches_actual_schedule_multiplier_and_is_radially_constant(reference) -> None:
    chi = actual_schedule_chi_coefficient_state(reference)
    for eta in (-0.31, -0.011, 0.27):
        assert chi.jet(0, 0, eta) == pytest.approx(reference.reference.chi0(eta), rel=2e-14, abs=2e-14)
        assert chi.jet(0, 1, eta) == pytest.approx(_chi_eta(reference, eta), rel=2e-9, abs=2e-12)
        for n in (1, 2, 5):
            assert chi.jet(n, 0, eta) == 0.0
            assert chi.jet(n, 2, eta) == 0.0


def test_natural_operator_is_exact_half_j2_chi_product_on_actual_reference(reference, natural_operator) -> None:
    eta = -0.17
    chi0 = reference.reference.chi0(eta)
    chi1 = _chi_eta(reference, eta)
    output = natural_operator(reference.phi)

    assert output.jet(0, 0, eta) == 0.0
    assert output.jet(0, 1, eta) == 0.0
    for n in range(1, 6):
        divisor = float(n * (n + 1))
        a0 = reference.phi.jet(n - 1, 0, eta)
        a1 = reference.phi.jet(n - 1, 1, eta)
        expected0 = 0.5 * chi0 * a0 / divisor
        expected1 = 0.5 * (chi1 * a0 + chi0 * a1) / divisor
        assert output.jet(n, 0, eta) == pytest.approx(expected0, rel=3e-13, abs=3e-15)
        assert output.jet(n, 1, eta) == pytest.approx(expected1, rel=3e-11, abs=3e-13)


def test_natural_operator_matches_reference_resolvent_recurrence(reference, natural_operator) -> None:
    """The landed closed-form reference must satisfy phi0 + Q(phi0) = one."""

    eta = 0.23
    q_phi = natural_operator(reference.phi)
    assert reference.phi.jet(0, 0, eta) + q_phi.jet(0, 0, eta) == pytest.approx(1.0)
    for n in range(1, 7):
        assert reference.phi.jet(n, 0, eta) + q_phi.jet(n, 0, eta) == pytest.approx(
            0.0,
            rel=2e-12,
            abs=2e-15,
        )


def test_natural_operator_fails_closed_on_wrong_scale(reference, natural_operator) -> None:
    wrong = AxisCoefficientJetState(
        epsilon=2.0 * reference.epsilon,
        origin="wrong-scale regression state",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        natural_operator(wrong)


def test_natural_operator_truth_boundary_is_explicit(reference, natural_operator) -> None:
    assert natural_operator.epsilon == reference.epsilon
    assert natural_operator.paper_exact is False
    assert natural_operator.global_axis_norm_certified is False
    assert natural_operator.natural_resolvent_materialized is False


def test_natural_operator_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_natural_operator.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "naturalOperator" in layer["capability"]
    assert "naturalResolvent" in layer["remaining_boundary"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
