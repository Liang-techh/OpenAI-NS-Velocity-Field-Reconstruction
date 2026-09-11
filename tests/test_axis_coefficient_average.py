import json
from pathlib import Path

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.axis_coefficient_average import (
    average_jet,
    axis_coefficient_average,
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


def test_average_coefficients_match_independent_radial_integral(state) -> None:
    averaged = axis_coefficient_average(state.phi)
    eta = 0.31
    y = 0.43
    degree = 6

    # Independent oracle: integrate a finite actual-reference radial polynomial
    # in physical radial scale t.  Production acts coefficient-by-coefficient.
    def radial_polynomial(t: float) -> float:
        return sum(
            state.phi.coefficient(n, eta) * (t * y) ** n
            for n in range(degree + 1)
        )

    oracle, error = quad(radial_polynomial, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12)
    reconstructed = sum(
        averaged.coefficient(n, eta) * y**n for n in range(degree + 1)
    )
    assert error < 1e-9
    assert reconstructed == pytest.approx(oracle, rel=2e-11, abs=2e-11)


def test_average_parameter_jet_matches_finite_difference(state) -> None:
    averaged = axis_coefficient_average(state.phi)
    eta = -0.19
    n = 2
    step = 2.0e-5

    def direct_average_coefficient(x: float) -> float:
        return state.phi.coefficient(n, x) / float(n + 1)

    oracle = (
        direct_average_coefficient(eta + step)
        - direct_average_coefficient(eta - step)
    ) / (2.0 * step)
    assert averaged.jet(n, 1, eta) == pytest.approx(oracle, rel=2e-5, abs=2e-7)


def test_axial_reference_average_has_exact_half_degree_one(state) -> None:
    averaged = axis_coefficient_average(state.u)
    eta = 0.22
    for m in range(4):
        assert averaged.jet(1, m, eta) == pytest.approx(
            0.5 * state.u.jet(1, m, eta), rel=2e-14, abs=1e-12
        )
        assert averaged.jet(0, m, eta) == 0.0
        assert averaged.jet(2, m, eta) == 0.0


def test_average_jet_matches_pinned_row_formula(state) -> None:
    eta = -0.47
    for n in range(6):
        for m in range(3):
            assert average_jet(state.phi, n, m, eta) == pytest.approx(
                state.phi.jet(n, m, eta) / float(n + 1),
                rel=2e-15,
                abs=1e-14,
            )


def test_average_stays_fail_closed_for_axis_space_membership(state) -> None:
    averaged = axis_coefficient_average(state.phi)
    assert averaged.paper_exact is False
    assert averaged.global_axis_norm_certified is False
    assert averaged.epsilon == state.epsilon
    assert "AxisOperators average" in averaged.origin


def test_average_preserves_underlying_index_and_window_guards(state) -> None:
    averaged = axis_coefficient_average(state.phi)
    with pytest.raises(ValueError, match="nonnegative integer"):
        averaged.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        averaged.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        averaged.jet(0, 0, 1.100001)


def test_average_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_average.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.averageData" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
