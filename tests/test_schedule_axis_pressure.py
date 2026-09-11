import math

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.outgoing_tail import (
    OutgoingCoreParameters,
    TailData,
    final_angular,
)
from openai_ns_reconstruction.schedule_axis_pressure import (
    axis_pressure,
    axis_pressure_derivative,
    exponent_weighted_mass,
    ideal_prefix_mass,
    pressure_kernel,
    schedule_pressure_mass,
)


def _data() -> TailData:
    return TailData(OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0), h=0.01)


def _direct_mass(data: TailData, eta: float) -> float:
    """Independent SciPy integral of the unfactorized ``final_angular**2``."""

    # Splitting at the known smooth-transition boundaries helps the adaptive
    # oracle without copying the production analytic plateau reductions.
    points = (
        0.0,
        1.0,
        data.core.drop_length + 1.0,
        data.core.hold_start,
        data.core.endpoint,
        data.flatten_end,
        data.release_start,
        data.release_start + 1.0,
        data.release_start + data.second_ramp_start,
        data.release_start + data.ramp_end,
        data.tail_start + 1.0,
        data.tail_end,
    )
    fn = lambda y: final_angular(data, float(y), eta) ** 2
    total = quad(fn, -math.inf, points[0], epsabs=2e-11, epsrel=2e-11, limit=200)[0]
    for a, b in zip(points, points[1:]):
        total += quad(fn, a, b, epsabs=2e-11, epsrel=2e-11, limit=300)[0]
    total += quad(fn, points[-1], math.inf, epsabs=2e-11, epsrel=2e-11, limit=200)[0]
    return total


def test_pressure_kernel_is_stable_on_finite_real_axis() -> None:
    assert pressure_kernel(1.0, 0.0) == 1.0
    assert pressure_kernel(1.0, 2.0) == pytest.approx(1.0 / 25.0)
    assert pressure_kernel(0.0, 1e308) == 1.0
    value = pressure_kernel(0.4, 1e308)
    assert math.isfinite(value) and value >= 0.0
    with pytest.raises(ValueError):
        pressure_kernel(-0.01, 0.0)
    with pytest.raises(ValueError):
        pressure_kernel(1.01, 0.0)
    with pytest.raises(ValueError):
        pressure_kernel(0.5, math.inf)


def test_actual_schedule_pressure_matches_independent_improper_integral() -> None:
    data = _data()
    for eta in (0.0, 0.35, 1.0):
        direct = _direct_mass(data, eta)
        production = schedule_pressure_mass(data, eta)
        assert production == pytest.approx(direct, rel=3e-9, abs=2e-10)
        assert axis_pressure(data, eta) == pytest.approx(-0.5 * direct, rel=3e-9, abs=1e-10)


def test_ideal_prefix_and_lean_pressure_bound() -> None:
    data = _data()
    for eta in (0.0, 0.4, 1.0, 2.0):
        expected_prefix = 5.0 * data.core.P**2 / (1.0 + eta * eta) ** 2
        assert ideal_prefix_mass(data, eta) == pytest.approx(expected_prefix, rel=3e-15)
        # SchedulePressure.axisPressure_lower_bound: the positive-y schedule
        # only adds nonnegative pressure mass to the exact ideal prefix.
        lean_upper_bound = -2.5 * data.core.P**2 / (1.0 + eta * eta) ** 2
        assert axis_pressure(data, eta) <= lean_upper_bound


def test_evenness_and_derivative_formula_against_finite_difference() -> None:
    data = _data()
    for eta in (0.2, 0.7):
        assert axis_pressure(data, -eta) == pytest.approx(axis_pressure(data, eta), rel=2e-13)
        assert axis_pressure_derivative(data, -eta) == pytest.approx(
            -axis_pressure_derivative(data, eta), rel=2e-13
        )
        step = 2e-5
        finite_difference = (
            axis_pressure(data, eta + step) - axis_pressure(data, eta - step)
        ) / (2.0 * step)
        assert axis_pressure_derivative(data, eta) == pytest.approx(
            finite_difference, rel=2e-7, abs=2e-9
        )
    assert axis_pressure_derivative(data, 0.0) == 0.0
    assert exponent_weighted_mass(data, 0.0) > 0.0


def test_schedule_pressure_fails_closed_on_invalid_inputs() -> None:
    data = _data()
    with pytest.raises(TypeError):
        schedule_pressure_mass(object(), 0.0)
    with pytest.raises(ValueError):
        axis_pressure(data, math.nan)
    with pytest.raises(ValueError):
        axis_pressure_derivative(data, math.inf)
