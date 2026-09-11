import json
from pathlib import Path

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.axis_coefficient_primitive import (
    axis_coefficient_primitive,
    primitive_jet,
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


def test_primitive_coefficients_match_independent_radial_integral(state) -> None:
    primitive = axis_coefficient_primitive(state.phi)
    eta = 0.27
    y = 0.41
    degree = 6

    # Independent oracle: integrate the finite actual-reference radial
    # polynomial directly in the physical radial variable s.  Production uses
    # only the pinned coefficient row map.
    def radial_polynomial(s: float) -> float:
        return sum(
            state.phi.coefficient(n, eta) * s**n for n in range(degree + 1)
        )

    oracle, error = quad(radial_polynomial, 0.0, y, epsabs=1e-12, epsrel=1e-12)
    reconstructed = sum(
        primitive.coefficient(n, eta) * y**n for n in range(degree + 2)
    )
    assert error < 1e-9
    assert reconstructed == pytest.approx(oracle, rel=2e-11, abs=2e-11)


def test_primitive_parameter_jet_matches_finite_difference(state) -> None:
    primitive = axis_coefficient_primitive(state.phi)
    eta = -0.23
    n = 3
    step = 2.0e-5

    def direct_primitive_coefficient(x: float) -> float:
        return state.phi.coefficient(n - 1, x) / float(n)

    oracle = (
        direct_primitive_coefficient(eta + step)
        - direct_primitive_coefficient(eta - step)
    ) / (2.0 * step)
    assert primitive.jet(n, 1, eta) == pytest.approx(oracle, rel=2e-5, abs=2e-7)


def test_axial_reference_primitive_moves_degree_one_to_degree_two(state) -> None:
    primitive = axis_coefficient_primitive(state.u)
    eta = 0.18
    for m in range(4):
        assert primitive.jet(0, m, eta) == 0.0
        assert primitive.jet(1, m, eta) == 0.0
        assert primitive.jet(2, m, eta) == pytest.approx(
            0.5 * state.u.jet(1, m, eta), rel=2e-14, abs=1e-12
        )
        assert primitive.jet(3, m, eta) == 0.0


def test_primitive_jet_matches_pinned_row_formula(state) -> None:
    eta = -0.39
    for m in range(3):
        assert primitive_jet(state.phi, 0, m, eta) == 0.0
        for n in range(1, 7):
            assert primitive_jet(state.phi, n, m, eta) == pytest.approx(
                state.phi.jet(n - 1, m, eta) / float(n),
                rel=2e-15,
                abs=1e-14,
            )


def test_primitive_stays_fail_closed_for_axis_space_membership(state) -> None:
    primitive = axis_coefficient_primitive(state.phi)
    assert primitive.paper_exact is False
    assert primitive.global_axis_norm_certified is False
    assert primitive.epsilon == state.epsilon
    assert "AxisOperators primitive" in primitive.origin


def test_primitive_preserves_underlying_index_and_window_guards(state) -> None:
    primitive = axis_coefficient_primitive(state.phi)
    with pytest.raises(ValueError, match="nonnegative integer"):
        primitive.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        primitive.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        primitive.jet(0, 0, 1.100001)


def test_primitive_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_primitive.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.primitiveData" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
