import math

import pytest

from openai_ns_reconstruction.axis_reference_axial_jets import (
    actual_schedule_axial_reference_jets,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_axis_pressure import axis_pressure, axis_pressure_derivative
from openai_ns_reconstruction.schedule_axis_pressure_jets import (
    axis_pressure_normalized_taylor,
    pressure_kernel_normalized_taylor,
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def jets():
    return actual_schedule_axial_reference_jets(_schedule_data(), 0.05)


def test_kernel_taylor_has_independent_closed_form_specializations() -> None:
    # a=1: (1+t^2)^-2 = 1 - 2 t^2 + 3 t^4 + ... at eta=0.
    assert pressure_kernel_normalized_taylor(1.0, 0.0, 4) == pytest.approx(
        (1.0, 0.0, -2.0, 0.0, 3.0), rel=2e-15, abs=2e-15
    )
    # a=1/2: (1+t^2)^-1 = 1 - t^2 + t^4 + ...
    assert pressure_kernel_normalized_taylor(0.5, 0.0, 4) == pytest.approx(
        (1.0, 0.0, -1.0, 0.0, 1.0), rel=2e-15, abs=2e-15
    )


def test_pressure_taylor_matches_landed_value_and_first_derivative() -> None:
    data = _schedule_data()
    for eta in (-0.55, 0.0, 0.43):
        pressure_jet = axis_pressure_normalized_taylor(data, eta, 1)
        assert pressure_jet[0] == pytest.approx(axis_pressure(data, eta), rel=3e-11, abs=2e-10)
        assert pressure_jet[1] == pytest.approx(
            axis_pressure_derivative(data, eta), rel=3e-10, abs=2e-10
        )


def test_axial_zeroth_parameter_jet_matches_landed_reference_pair(jets) -> None:
    for eta in (-0.7, 0.0, 0.52):
        assert jets.parameter_jet(1, 0, eta) == pytest.approx(
            jets.reference.u_coefficient(1, eta), rel=3e-10, abs=2e-10
        )
        for n in (0, 2, 4):
            assert jets.parameter_jet(n, 0, eta) == 0.0


def test_first_axial_parameter_jet_matches_independent_finite_difference(jets) -> None:
    eta = 0.31
    step = 2e-5
    finite_difference = (
        jets.reference.u_coefficient(1, eta + step)
        - jets.reference.u_coefficient(1, eta - step)
    ) / (2.0 * step)
    assert jets.parameter_jet(1, 1, eta) == pytest.approx(
        finite_difference, rel=8e-7, abs=2e-7
    )


def test_second_axial_parameter_jet_matches_independent_second_difference(jets) -> None:
    eta = -0.27
    step = 2e-4
    center = jets.reference.u_coefficient(1, eta)
    second_difference = (
        jets.reference.u_coefficient(1, eta + step)
        - 2.0 * center
        + jets.reference.u_coefficient(1, eta - step)
    ) / (step * step)
    assert jets.parameter_jet(1, 2, eta) == pytest.approx(
        second_difference, rel=2e-4, abs=5e-5
    )


def test_axial_jets_use_actual_schedule_epsilon_and_fail_closed(jets) -> None:
    assert jets.epsilon == jets.analytic_inputs.neighborhood.radius / 2.0
    assert jets.epsilon > 0.0
    assert jets.paper_exact is False

    n, m, eta = 1, 1, 0.17
    coordinate = jets.normalized_axis_coordinate(n, m, eta)
    assert math.isfinite(coordinate)

    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(1, -1, 0.0)
    with pytest.raises(ValueError, match="pinned coefficient window"):
        jets.parameter_jet(1, 0, 1.100001)
