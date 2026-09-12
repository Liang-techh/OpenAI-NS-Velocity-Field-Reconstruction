import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_inverse_dot_product import (
    axis_coefficient_inverse_dot_product,
    inverse_dot_product_jet,
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


def _radial_dot_view(source: AxisCoefficientJetState) -> AxisCoefficientJetState:
    """Exact test-side view of ``Y * partial_Y source`` on radial coefficients."""

    return AxisCoefficientJetState(
        epsilon=source.epsilon,
        origin=f"exact radial-dot view of ({source.origin})",
        _jet_provider=lambda n, m, eta: float(n) * source.jet(n, m, eta),
    )


@pytest.mark.parametrize("r", [1, 2])
def test_inverse_dot_product_matches_independent_operator_composition(state, r: int) -> None:
    direct = axis_coefficient_inverse_dot_product(state.phi, state.u, r)
    dotted = _radial_dot_view(state.u)
    composed = axis_coefficient_regular_inverse(
        axis_coefficient_product(state.phi, dotted),
        r,
    )
    eta = -0.24

    # Production uses the direct pinned differentialFamily specialization.
    # The oracle travels through two separately landed/tested primitives:
    # productFamily followed by regularInverse, with a test-side exact Euler
    # radial derivative on the second operand.
    for n in range(8):
        for m in range(3):
            assert direct.jet(n, m, eta) == pytest.approx(
                composed.jet(n, m, eta),
                rel=3e-14,
                abs=2e-11,
            )


@pytest.mark.parametrize("r", [1, 2])
def test_actual_axial_square_has_exact_low_order_inverse_dot_identity(state, r: int) -> None:
    result = axis_coefficient_inverse_dot_product(state.u, state.u, r)
    eta = 0.19
    divisor = radial_divisor(r, 2)
    u0 = state.u.jet(1, 0, eta)
    u1 = state.u.jet(1, 1, eta)
    u2 = state.u.jet(1, 2, eta)

    # The actual axial reference has only radial degree one.  Since
    # Y*partial_Y acts as multiplication by one on that row, the source is
    # u^2 at radial degree two and its regular inverse lives at degree three.
    assert result.jet(3, 0, eta) == pytest.approx(u0 * u0 / divisor, rel=3e-14, abs=1e-11)
    assert result.jet(3, 1, eta) == pytest.approx(
        2.0 * u0 * u1 / divisor,
        rel=3e-13,
        abs=1e-8,
    )
    assert result.jet(3, 2, eta) == pytest.approx(
        (2.0 * u0 * u2 + 2.0 * u1 * u1) / divisor,
        rel=5e-12,
        abs=2e-6,
    )
    for n in (0, 1, 2, 4, 5):
        assert result.jet(n, 0, eta) == pytest.approx(0.0, abs=1e-13)


def test_inverse_dot_product_parameter_jet_matches_finite_difference(state) -> None:
    eta = -0.16
    step = 2.0e-5
    n = 4
    r = 2

    def direct_value(x: float) -> float:
        return inverse_dot_product_jet(state.phi, state.u, r, n, 0, x)

    oracle = (direct_value(eta + step) - direct_value(eta - step)) / (2.0 * step)
    actual = inverse_dot_product_jet(state.phi, state.u, r, n, 1, eta)
    assert actual == pytest.approx(oracle, rel=3e-5, abs=2e-7)


def test_inverse_dot_product_satisfies_truncated_physical_radial_equation(state) -> None:
    result = axis_coefficient_inverse_dot_product(state.phi, state.u, 2)
    eta = 0.14
    y = 0.31
    degree = 6

    # Independent source oracle: reconstruct the zeroth-eta-jet polynomials
    # and form phi(Y) * (Y * partial_Y u(Y)) directly.  Production instead
    # evaluates the pinned differentialFamily double sum.
    source_coefficients = []
    for q in range(degree + 1):
        coefficient = 0.0
        for i in range(q + 1):
            j = q - i
            coefficient += (
                state.phi.coefficient(i, eta)
                * float(j)
                * state.u.coefficient(j, eta)
            )
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
        rel=2e-12,
        abs=5e-10,
    )


def test_inverse_dot_product_stays_fail_closed_for_axis_space_membership(state) -> None:
    result = axis_coefficient_inverse_dot_product(state.phi, state.u, 1)
    assert result.paper_exact is False
    assert result.global_axis_norm_certified is False
    assert result.epsilon == state.epsilon
    assert "AxisOperators inverseDotProduct(r=1)" in result.origin


def test_inverse_dot_product_rejects_invalid_inputs(state) -> None:
    other = AxisCoefficientJetState(
        epsilon=2.0 * state.epsilon,
        origin="mismatched test state",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        axis_coefficient_inverse_dot_product(state.phi, other, 1)

    for bad_r in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            axis_coefficient_inverse_dot_product(state.phi, state.u, bad_r)

    result = axis_coefficient_inverse_dot_product(state.phi, state.u, 2)
    with pytest.raises(ValueError, match="nonnegative integer"):
        result.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        result.jet(True, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        result.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        result.jet(0, 0, 1.100001)


def test_inverse_dot_product_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_inverse_dot_product.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.inverseDotProduct" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
