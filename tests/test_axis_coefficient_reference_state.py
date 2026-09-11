import json
import math
from pathlib import Path

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.axis_coefficient_reference_state import (
    WINDOW_LEFT,
    WINDOW_RIGHT,
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


def test_reference_state_uses_one_actual_schedule_scale_and_stays_fail_closed(state) -> None:
    assert state.epsilon == state.angular_jets.epsilon == state.axial_jets.epsilon
    assert state.epsilon > 0.0
    assert state.reference.sigma == state.angular_jets.reference.sigma
    assert state.reference.sigma == state.axial_jets.reference.sigma
    assert state.paper_exact is False
    assert state.phi.paper_exact is False
    assert state.u.paper_exact is False
    assert state.global_axis_norm_certified is False
    assert state.phi.global_axis_norm_certified is False
    assert state.u.global_axis_norm_certified is False


def test_reference_state_zeroth_jets_match_landed_reference_pair(state) -> None:
    for eta in (-0.83, -0.17, 0.0, 0.61):
        for n in range(5):
            phi, u = state.coefficient_pair(n, eta)
            assert phi == pytest.approx(
                state.reference.phi_coefficient(n, eta), rel=3e-13, abs=1e-300
            )
            assert u == pytest.approx(
                state.reference.u_coefficient(n, eta), rel=3e-10, abs=2e-10
            )


def test_reference_state_exposes_actual_derivatives_not_independent_arrays(state) -> None:
    eta = 0.23
    for n, m in ((0, 0), (1, 1), (3, 2)):
        assert state.phi.jet(n, m, eta) == pytest.approx(
            state.angular_jets.parameter_jet(n, m, eta), rel=2e-14, abs=1e-300
        )
    for n, m in ((0, 2), (1, 0), (1, 2), (4, 1)):
        assert state.u.jet(n, m, eta) == pytest.approx(
            state.axial_jets.parameter_jet(n, m, eta), rel=3e-12, abs=2e-10
        )


def test_both_reference_components_satisfy_independent_ftc_checks(state) -> None:
    left = -0.21
    right = 0.26

    angular_integral, _ = quad(
        lambda eta: state.phi.jet(2, 2, eta),
        left,
        right,
        epsabs=2e-11,
        epsrel=2e-11,
        limit=100,
    )
    angular_difference = state.phi.jet(2, 1, right) - state.phi.jet(2, 1, left)
    assert angular_integral == pytest.approx(angular_difference, rel=2e-9, abs=2e-11)

    # The axial path is deliberately checked through an independent adaptive
    # quadrature of the first actual eta derivative.  Production obtains these
    # derivatives by truncated Taylor algebra through SchedulePressure.
    axial_integral, _ = quad(
        lambda eta: state.u.jet(1, 1, eta),
        left,
        right,
        epsabs=2e-8,
        epsrel=2e-8,
        limit=80,
    )
    axial_difference = state.u.jet(1, 0, right) - state.u.jet(1, 0, left)
    assert axial_integral == pytest.approx(axial_difference, rel=2e-6, abs=3e-8)


def test_normalized_coordinate_matches_pinned_weight_formula(state) -> None:
    n, m, eta = 2, 1, 0.19
    weight = (
        (1.0 / 20.0) ** n
        * state.epsilon ** (-m)
        * math.factorial(m)
        * math.comb(n + m, m)
        / ((n + 1.0) ** 2 * (m + 1.0) ** 2)
    )
    assert state.phi.normalized_coordinate(n, m, eta) == pytest.approx(
        state.phi.jet(n, m, eta) / weight, rel=2e-15
    )


def test_reference_state_uses_exact_pinned_window_guards(state) -> None:
    assert WINDOW_LEFT == -1.1
    assert WINDOW_RIGHT == 1.1
    state.phi.jet(0, 0, WINDOW_LEFT)
    state.u.jet(1, 0, WINDOW_RIGHT)

    with pytest.raises(ValueError, match="pinned window"):
        state.phi.jet(0, 0, WINDOW_RIGHT + 1e-6)
    with pytest.raises(ValueError, match="pinned window"):
        state.u.jet(1, 0, WINDOW_LEFT - 1e-6)
    with pytest.raises(ValueError, match="nonnegative integer"):
        state.phi.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        state.u.jet(1, -1, 0.0)


def test_reference_state_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_reference_state.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
