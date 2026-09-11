import numpy as np
import pytest

from openai_ns_reconstruction.past_extension import (
    past_pressure,
    past_residual_numeric,
    past_velocity,
    zero_before_scalar,
    zero_before_vector,
)
from openai_ns_reconstruction.time_localization import (
    activated_residual_formula_numeric,
)


def base_velocity(x: float, y: float, z: float, t: float) -> np.ndarray:
    # Smooth divergence-free solid rotation with time-dependent amplitude.
    return (1.0 + t) * np.array([-y, x, 0.0])


def base_pressure(x: float, y: float, z: float, t: float) -> float:
    return (1.0 - 0.25 * t) * (x + 2.0 * y - 0.5 * z)


def test_zero_before_short_circuits_negative_time() -> None:
    def unavailable_vector(x: float, y: float, z: float, t: float) -> np.ndarray:
        raise AssertionError("negative branch must not sample the input vector field")

    def unavailable_scalar(x: float, y: float, z: float, t: float) -> float:
        raise AssertionError("negative branch must not sample the input scalar field")

    assert np.array_equal(
        zero_before_vector(unavailable_vector)(0.2, -0.1, 0.3, -4.0),
        np.zeros(3),
    )
    assert zero_before_scalar(unavailable_scalar)(0.2, -0.1, 0.3, -4.0) == 0.0


def test_past_fields_have_the_pinned_zero_germ_and_late_identity() -> None:
    v = past_velocity(base_velocity)
    q = past_pressure(base_pressure)

    for t in (-2.0, -0.1, 0.0, 0.2, 3.0 / 8.0):
        assert np.array_equal(v(0.3, -0.4, 0.2, t), np.zeros(3))
        assert q(0.3, -0.4, 0.2, t) == 0.0

    point = (0.3, -0.4, 0.2, 0.8)
    assert np.array_equal(v(*point), base_velocity(*point))
    assert q(*point) == base_pressure(*point)


def test_past_residual_is_exactly_zero_nonpositive_without_sampling_inputs() -> None:
    def unavailable_vector(x: float, y: float, z: float, t: float) -> np.ndarray:
        raise AssertionError("nonpositive residual branch must not sample velocity")

    def unavailable_scalar(x: float, y: float, z: float, t: float) -> float:
        raise AssertionError("nonpositive residual branch must not sample pressure")

    for t in (-10.0, -0.2, 0.0):
        assert np.array_equal(
            past_residual_numeric(unavailable_vector, unavailable_scalar, 0.1, 0.2, -0.3, t),
            np.zeros(3),
        )


def test_past_residual_matches_independent_activation_formula_on_transition() -> None:
    point = (0.31, -0.27, 0.19, 0.55)
    direct = past_residual_numeric(
        base_velocity,
        base_pressure,
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


def test_past_residual_respects_open_singular_endpoint() -> None:
    with pytest.raises(ValueError, match="t<endpoint"):
        past_residual_numeric(base_velocity, base_pressure, 0.1, 0.2, 0.3, 1.0)
    with pytest.raises(ValueError, match="positive"):
        past_residual_numeric(
            base_velocity, base_pressure, 0.1, 0.2, 0.3, 0.5, endpoint=0.0
        )
