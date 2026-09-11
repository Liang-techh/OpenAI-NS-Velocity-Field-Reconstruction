import json
from pathlib import Path

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.axis_coefficient_parameter_primitive import (
    axis_coefficient_parameter_primitive,
    parameter_primitive_jet,
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


def test_parameter_primitive_matches_independent_eta_difference_and_radial_integral(state) -> None:
    parameter_primitive = axis_coefficient_parameter_primitive(state.phi)
    eta = 0.19
    y = 0.37
    degree = 6
    step = 1.5e-5

    # Independent oracle: differentiate the actual-reference radial polynomial
    # by a centered eta difference, then integrate that derivative in the
    # physical radial variable.  Production instead reads the analytic m+1 jet
    # and applies the pinned row map.
    def radial_polynomial(s: float, x: float) -> float:
        return sum(
            state.phi.coefficient(n, x) * s**n for n in range(degree + 1)
        )

    def eta_difference(s: float) -> float:
        return (
            radial_polynomial(s, eta + step)
            - radial_polynomial(s, eta - step)
        ) / (2.0 * step)

    oracle, error = quad(eta_difference, 0.0, y, epsabs=1e-10, epsrel=1e-10)
    reconstructed = sum(
        parameter_primitive.coefficient(n, eta) * y**n
        for n in range(degree + 2)
    )
    assert error < 1e-7
    assert reconstructed == pytest.approx(oracle, rel=3e-5, abs=3e-7)


def test_parameter_primitive_higher_parameter_jet_matches_finite_difference(state) -> None:
    parameter_primitive = axis_coefficient_parameter_primitive(state.phi)
    eta = -0.21
    n = 4
    step = 2.0e-5

    def direct_parameter_primitive_coefficient(x: float) -> float:
        return state.phi.jet(n - 1, 1, x) / float(n)

    oracle = (
        direct_parameter_primitive_coefficient(eta + step)
        - direct_parameter_primitive_coefficient(eta - step)
    ) / (2.0 * step)
    assert parameter_primitive.jet(n, 1, eta) == pytest.approx(
        oracle, rel=3e-5, abs=3e-7
    )


def test_axial_reference_parameter_primitive_uses_next_eta_jet(state) -> None:
    parameter_primitive = axis_coefficient_parameter_primitive(state.u)
    eta = 0.18
    for m in range(4):
        assert parameter_primitive.jet(0, m, eta) == 0.0
        assert parameter_primitive.jet(1, m, eta) == 0.0
        assert parameter_primitive.jet(2, m, eta) == pytest.approx(
            0.5 * state.u.jet(1, m + 1, eta), rel=2e-14, abs=1e-12
        )
        assert parameter_primitive.jet(3, m, eta) == 0.0


def test_parameter_primitive_jet_matches_pinned_shift_formula(state) -> None:
    eta = -0.39
    for m in range(3):
        assert parameter_primitive_jet(state.phi, 0, m, eta) == 0.0
        for n in range(1, 7):
            assert parameter_primitive_jet(state.phi, n, m, eta) == pytest.approx(
                state.phi.jet(n - 1, m + 1, eta) / float(n),
                rel=2e-15,
                abs=1e-14,
            )


def test_parameter_primitive_stays_fail_closed_for_axis_space_membership(state) -> None:
    parameter_primitive = axis_coefficient_parameter_primitive(state.phi)
    assert parameter_primitive.paper_exact is False
    assert parameter_primitive.global_axis_norm_certified is False
    assert parameter_primitive.epsilon == state.epsilon
    assert "AxisOperators parameterPrimitive" in parameter_primitive.origin


def test_parameter_primitive_preserves_underlying_index_and_window_guards(state) -> None:
    parameter_primitive = axis_coefficient_parameter_primitive(state.phi)
    with pytest.raises(ValueError, match="nonnegative integer"):
        parameter_primitive.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        parameter_primitive.jet(True, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        parameter_primitive.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        parameter_primitive.jet(0, 0, 1.100001)


def test_parameter_primitive_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_parameter_primitive.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.parameterPrimitive" in layer["capability"]
    assert "m+1" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
