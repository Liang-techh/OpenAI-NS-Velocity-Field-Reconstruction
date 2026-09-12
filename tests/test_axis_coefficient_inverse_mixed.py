import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_inverse_mixed import (
    axis_coefficient_inverse_mixed,
    inverse_mixed_jet,
)
from openai_ns_reconstruction.axis_coefficient_product import axis_coefficient_product
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.axis_coefficient_regular_inverse import (
    axis_coefficient_regular_inverse,
    radial_divisor,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


ROOT = Path(__file__).resolve().parents[1]


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def state():
    return actual_schedule_reference_axis_state(_schedule_data(), 0.05)


def _parameter_derivative_view(source: AxisCoefficientJetState) -> AxisCoefficientJetState:
    """Exact test-side view of ``partial_eta source`` using landed jets."""

    return AxisCoefficientJetState(
        epsilon=source.epsilon,
        origin=f"exact parameter derivative view of ({source.origin})",
        _jet_provider=lambda n, m, eta: source.jet(n, m + 1, eta),
    )


def _radial_dot_view(source: AxisCoefficientJetState) -> AxisCoefficientJetState:
    """Exact test-side view of ``Y * partial_Y source`` on radial coefficients."""

    return AxisCoefficientJetState(
        epsilon=source.epsilon,
        origin=f"exact radial-dot view of ({source.origin})",
        _jet_provider=lambda n, m, eta: float(n) * source.jet(n, m, eta),
    )


@pytest.mark.parametrize("r", [1, 2])
def test_inverse_mixed_matches_independent_operator_composition(state, r: int) -> None:
    direct = axis_coefficient_inverse_mixed(state.u, state.phi, r)
    derivative = _parameter_derivative_view(state.u)
    dotted = _radial_dot_view(state.phi)
    composed = axis_coefficient_regular_inverse(
        axis_coefficient_product(derivative, dotted),
        r,
    )
    eta = -0.24

    # Production evaluates the direct pinned differentialFamily specialization.
    # The oracle travels through two separately landed/tested primitives after
    # exact test-side parameter/radial transformations of the actual states.
    for n in range(8):
        for m in range(3):
            assert direct.jet(n, m, eta) == pytest.approx(
                composed.jet(n, m, eta),
                rel=4e-14,
                abs=3e-11,
            )


@pytest.mark.parametrize("r", [1, 2])
def test_actual_axial_square_has_exact_low_order_inverse_mixed_identity(state, r: int) -> None:
    result = axis_coefficient_inverse_mixed(state.u, state.u, r)
    eta = 0.19
    divisor = radial_divisor(r, 2)
    u0 = state.u.jet(1, 0, eta)
    u1 = state.u.jet(1, 1, eta)
    u2 = state.u.jet(1, 2, eta)

    # The actual axial reference has only radial degree one, so its Euler
    # radial derivative equals itself.  J_r((partial_eta u) D_Y u) therefore
    # has only radial degree three and admits these closed parameter jets.
    assert result.jet(3, 0, eta) == pytest.approx(u1 * u0 / divisor, rel=3e-14, abs=1e-11)
    assert result.jet(3, 1, eta) == pytest.approx(
        (u2 * u0 + u1 * u1) / divisor,
        rel=3e-13,
        abs=1e-8,
    )
    for n in (0, 1, 2, 4, 5):
        assert result.jet(n, 0, eta) == pytest.approx(0.0, abs=1e-13)


def test_inverse_mixed_parameter_jet_matches_finite_difference(state) -> None:
    eta = -0.16
    step = 2.0e-5
    n = 4
    r = 2

    def direct_value(x: float) -> float:
        return inverse_mixed_jet(state.u, state.phi, r, n, 0, x)

    oracle = (direct_value(eta + step) - direct_value(eta - step)) / (2.0 * step)
    actual = inverse_mixed_jet(state.u, state.phi, r, n, 1, eta)
    assert actual == pytest.approx(oracle, rel=4e-5, abs=3e-7)


def test_inverse_mixed_satisfies_truncated_physical_radial_equation(state) -> None:
    result = axis_coefficient_inverse_mixed(state.u, state.phi, 2)
    eta = 0.14
    step = 2.0e-5
    y = 0.31
    degree = 5

    # Independent source oracle: finite-difference the first factor in eta and
    # apply the Euler radial derivative to the second factor directly at the
    # coefficient level.  Production never finite-differences eta and instead
    # evaluates the pinned differentialFamily double sum.
    source_coefficients = []
    for q in range(degree + 1):
        coefficient = 0.0
        for i in range(q + 1):
            j = q - i
            left_eta = (
                state.u.coefficient(i, eta + step)
                - state.u.coefficient(i, eta - step)
            ) / (2.0 * step)
            coefficient += left_eta * float(j) * state.phi.coefficient(j, eta)
        source_coefficients.append(coefficient)

    source_value = sum(c * y**q for q, c in enumerate(source_coefficients))
    first = sum(
        n * result.coefficient(n, eta) * y ** (n - 1)
        for n in range(1, degree + 2)
    )
    second = sum(
        n * (n - 1) * result.coefficient(n, eta) * y ** (n - 2)
        for n in range(2, degree + 2)
    )
    assert y * second + 2.0 * first == pytest.approx(
        source_value,
        rel=5e-5,
        abs=8e-7,
    )


def test_inverse_mixed_stays_fail_closed_for_axis_space_membership(state) -> None:
    result = axis_coefficient_inverse_mixed(state.u, state.phi, 1)
    assert result.paper_exact is False
    assert result.global_axis_norm_certified is False
    assert result.epsilon == state.epsilon
    assert "AxisOperators inverseMixed(r=1)" in result.origin


def test_inverse_mixed_rejects_invalid_inputs(state) -> None:
    other = AxisCoefficientJetState(
        epsilon=2.0 * state.epsilon,
        origin="mismatched test state",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        axis_coefficient_inverse_mixed(state.u, other, 1)

    for bad_r in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            axis_coefficient_inverse_mixed(state.u, state.phi, bad_r)

    result = axis_coefficient_inverse_mixed(state.u, state.phi, 2)
    with pytest.raises(ValueError, match="nonnegative integer"):
        result.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        result.jet(True, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        result.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        result.jet(0, 0, 1.100001)


def test_inverse_mixed_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_inverse_mixed.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.inverseMixed" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert "coefficientOperators" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
