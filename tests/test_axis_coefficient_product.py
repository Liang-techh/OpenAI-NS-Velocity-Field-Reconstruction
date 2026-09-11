import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_product import (
    axis_coefficient_product,
    product_jet,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    actual_schedule_reference_axis_state,
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


def test_axial_square_matches_independent_low_order_leibniz_identities(state) -> None:
    product = axis_coefficient_product(state.u, state.u)
    eta = 0.23
    u0 = state.u.jet(1, 0, eta)
    u1 = state.u.jet(1, 1, eta)
    u2 = state.u.jet(1, 2, eta)

    assert product.jet(2, 0, eta) == pytest.approx(u0 * u0, rel=2e-14, abs=1e-12)
    assert product.jet(2, 1, eta) == pytest.approx(2.0 * u0 * u1, rel=2e-14, abs=1e-10)
    assert product.jet(2, 2, eta) == pytest.approx(
        2.0 * u1 * u1 + 2.0 * u0 * u2,
        rel=3e-13,
        abs=1e-8,
    )

    # The actual axial reference has only radial degree one.  Its square must
    # therefore have exactly radial degree two at every parameter-jet order.
    for n in (0, 1, 3, 4):
        assert product.jet(n, 2, eta) == 0.0


def test_mixed_product_parameter_jet_matches_finite_difference_of_values(state) -> None:
    product = axis_coefficient_product(state.phi, state.u)
    eta = -0.17
    radial_degree = 3
    step = 2.0e-5

    # Independent oracle: differentiate the zeroth-jet radial convolution by a
    # centered finite difference.  Production uses the exact Leibniz jet sum.
    def direct_value(x: float) -> float:
        return sum(
            state.phi.coefficient(i, x) * state.u.coefficient(radial_degree - i, x)
            for i in range(radial_degree + 1)
        )

    oracle = (direct_value(eta + step) - direct_value(eta - step)) / (2.0 * step)
    assert product.jet(radial_degree, 1, eta) == pytest.approx(
        oracle,
        rel=2e-5,
        abs=5e-7,
    )


def test_product_zeroth_jet_is_direct_radial_convolution(state) -> None:
    eta = 0.41
    for n in range(6):
        expected = sum(
            state.phi.coefficient(i, eta) * state.u.coefficient(n - i, eta)
            for i in range(n + 1)
        )
        assert product_jet(state.phi, state.u, n, 0, eta) == pytest.approx(
            expected,
            rel=3e-14,
            abs=1e-12,
        )


def test_product_stays_fail_closed_for_axis_space_membership(state) -> None:
    product = axis_coefficient_product(state.phi, state.u)
    assert product.paper_exact is False
    assert product.global_axis_norm_certified is False
    assert "AxisOperators product" in product.origin


def test_product_rejects_mixed_weight_scales(state) -> None:
    other = AxisCoefficientJetState(
        epsilon=2.0 * state.epsilon,
        origin="mismatched test state",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        axis_coefficient_product(state.phi, other)
    with pytest.raises(ValueError, match="matching epsilon"):
        product_jet(state.phi, other, 0, 0, 0.0)


def test_product_preserves_underlying_index_and_window_guards(state) -> None:
    product = axis_coefficient_product(state.phi, state.u)
    with pytest.raises(ValueError, match="nonnegative integer"):
        product.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        product.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        product.jet(0, 0, 1.100001)


def test_product_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_product.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.productFamily" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
