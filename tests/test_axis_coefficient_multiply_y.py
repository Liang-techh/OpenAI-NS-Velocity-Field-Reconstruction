import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_multiply_y import (
    axis_coefficient_multiply_y,
    multiply_y_jet,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
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


def test_multiply_y_matches_independent_polynomial_product(state) -> None:
    multiplied = axis_coefficient_multiply_y(state.phi)
    eta = 0.27
    y = 0.41
    degree = 6

    # Independent physical-polynomial oracle: multiply the finite truncation by
    # the radial variable y directly. Production uses only the pinned row map.
    source_value = sum(
        state.phi.coefficient(n, eta) * y**n for n in range(degree + 1)
    )
    oracle = y * source_value
    reconstructed = sum(
        multiplied.coefficient(n, eta) * y**n for n in range(degree + 2)
    )
    assert reconstructed == pytest.approx(oracle, rel=2e-13, abs=2e-13)


def test_multiply_y_parameter_jet_matches_finite_difference(state) -> None:
    multiplied = axis_coefficient_multiply_y(state.phi)
    eta = -0.23
    n = 3
    step = 2.0e-5

    def direct_multiplied_coefficient(x: float) -> float:
        return state.phi.coefficient(n - 1, x)

    oracle = (
        direct_multiplied_coefficient(eta + step)
        - direct_multiplied_coefficient(eta - step)
    ) / (2.0 * step)
    assert multiplied.jet(n, 1, eta) == pytest.approx(oracle, rel=2e-5, abs=2e-7)


def test_axial_reference_multiply_y_moves_degree_one_to_degree_two(state) -> None:
    multiplied = axis_coefficient_multiply_y(state.u)
    eta = 0.18
    for m in range(4):
        assert multiplied.jet(0, m, eta) == 0.0
        assert multiplied.jet(1, m, eta) == 0.0
        assert multiplied.jet(2, m, eta) == pytest.approx(
            state.u.jet(1, m, eta), rel=2e-14, abs=1e-12
        )
        assert multiplied.jet(3, m, eta) == 0.0


def test_multiply_y_jet_matches_pinned_row_formula(state) -> None:
    eta = -0.39
    for m in range(3):
        assert multiply_y_jet(state.phi, 0, m, eta) == 0.0
        for n in range(1, 7):
            assert multiply_y_jet(state.phi, n, m, eta) == pytest.approx(
                state.phi.jet(n - 1, m, eta), rel=2e-15, abs=1e-14
            )


def test_multiply_y_stays_fail_closed_for_axis_space_membership(state) -> None:
    multiplied = axis_coefficient_multiply_y(state.phi)
    assert multiplied.paper_exact is False
    assert multiplied.global_axis_norm_certified is False
    assert multiplied.epsilon == state.epsilon
    assert "AxisOperators mulY" in multiplied.origin


def test_multiply_y_preserves_underlying_index_and_window_guards(state) -> None:
    multiplied = axis_coefficient_multiply_y(state.phi)
    with pytest.raises(ValueError, match="nonnegative integer"):
        multiplied.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        multiplied.jet(True, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        multiplied.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        multiplied.jet(0, 0, 1.100001)


def test_multiply_y_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_multiply_y.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.multiplyYData" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
