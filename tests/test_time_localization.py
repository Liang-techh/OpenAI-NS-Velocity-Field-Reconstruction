import numpy as np

from openai_ns_reconstruction.time_localization import (
    EARLY_HALF_WIDTH,
    LATE_START,
    activated_pressure,
    activated_residual_formula_numeric,
    activated_velocity,
    time_switch,
    time_switch_derivative,
)
from openai_ns_reconstruction.verify import navier_stokes_residual_numeric


def base_velocity(x: float, y: float, z: float, t: float) -> np.ndarray:
    # Divergence-free, time-dependent solid rotation.
    return (1.0 + t) * np.array([-y, x, 0.0])


def base_pressure(x: float, y: float, z: float, t: float) -> float:
    return (1.0 - 0.25 * t) * (x + 2.0 * y - 0.5 * z)


def test_time_switch_has_official_early_and_late_geometry() -> None:
    for t in (-EARLY_HALF_WIDTH, -0.2, 0.0, 0.2, EARLY_HALF_WIDTH):
        assert time_switch(t) == 0.0
    for t in (LATE_START, 0.9, 1.0, 3.0):
        assert time_switch(t) == 1.0
    assert 0.0 < time_switch(0.55) < 1.0


def test_time_switch_derivative_matches_independent_finite_difference() -> None:
    t = 0.55
    h = 1e-6
    finite_difference = (time_switch(t + h) - time_switch(t - h)) / (2.0 * h)
    assert np.isclose(time_switch_derivative(t), finite_difference, rtol=2e-6, atol=2e-8)
    assert time_switch_derivative(t) > 0.0
    assert time_switch_derivative(0.0) == 0.0
    assert time_switch_derivative(0.9) == 0.0


def test_early_activation_short_circuits_unresolved_input_fields() -> None:
    def unavailable_velocity(x: float, y: float, z: float, t: float) -> np.ndarray:
        raise AssertionError("early switch should not evaluate the unresolved velocity")

    def unavailable_pressure(x: float, y: float, z: float, t: float) -> float:
        raise AssertionError("early switch should not evaluate the unresolved pressure")

    assert np.array_equal(activated_velocity(unavailable_velocity)(0.2, -0.1, 0.4, 0.0), np.zeros(3))
    assert activated_pressure(unavailable_pressure)(0.2, -0.1, 0.4, 0.0) == 0.0


def test_late_activation_is_exact_identity() -> None:
    point = (0.3, -0.4, 0.2, 0.8)
    assert np.array_equal(activated_velocity(base_velocity)(*point), base_velocity(*point))
    assert activated_pressure(base_pressure)(*point) == base_pressure(*point)


def test_activated_residual_formula_matches_independent_residual() -> None:
    point = (0.31, -0.27, 0.19, 0.55)
    v = activated_velocity(base_velocity)
    q = activated_pressure(base_pressure)
    direct = navier_stokes_residual_numeric(
        v,
        q,
        *point,
        viscosity=0.0,
        eps_space=2e-6,
        eps_time=2e-6,
    )
    formula = activated_residual_formula_numeric(
        base_velocity,
        base_pressure,
        *point,
        viscosity=0.0,
        eps_space=2e-6,
        eps_time=2e-6,
    )
    assert np.allclose(direct, formula, rtol=2e-5, atol=2e-7)
