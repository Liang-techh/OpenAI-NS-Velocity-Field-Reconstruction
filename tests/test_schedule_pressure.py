import math

import numpy as np
import pytest

from openai_ns_reconstruction.schedule_pressure import (
    EXPLICIT_FLATTEN_LENGTH,
    EXPLICIT_STEP_BOUND,
    flat_edge_unit,
    flatten_end,
    outgoing_sigma,
    outgoing_sigma_derivative,
    schedule_shape_exponent,
)


def _direct_edge_ratio(x: float) -> float:
    a = math.exp(-1.0 / (x * x))
    b = math.exp(-1.0 / ((1.0 - x) * (1.0 - x)))
    return a / (a + b)


def test_outgoing_sigma_matches_pinned_flat_edge_ratio() -> None:
    for x in (0.2, 0.25, 0.5, 0.75, 0.8):
        assert outgoing_sigma(x) == pytest.approx(_direct_edge_ratio(x), rel=2e-15, abs=0.0)
    assert outgoing_sigma(0.5) == pytest.approx(0.5, abs=2e-15)
    assert outgoing_sigma(-3.0) == 0.0
    assert outgoing_sigma(0.0) == 0.0
    assert outgoing_sigma(1.0) == 1.0
    assert outgoing_sigma(4.0) == 1.0


def test_flat_edge_is_the_squared_denominator_gaussian_edge() -> None:
    assert flat_edge_unit(-1.0) == 0.0
    assert flat_edge_unit(0.0) == 0.0
    assert flat_edge_unit(0.5) == pytest.approx(math.exp(-4.0), rel=2e-15)
    assert flat_edge_unit(2.0) == pytest.approx(math.exp(-0.25), rel=2e-15)


def test_outgoing_sigma_symmetry_monotonicity_and_derivative() -> None:
    xs = np.linspace(0.02, 0.98, 97)
    values = np.array([outgoing_sigma(float(x)) for x in xs])
    assert np.all(np.diff(values) >= -2e-15)
    for x in (0.1, 0.2, 0.37, 0.5, 0.71, 0.9):
        assert outgoing_sigma(x) + outgoing_sigma(1.0 - x) == pytest.approx(1.0, abs=3e-15)

    eps = 1e-6
    for x in (0.2, 0.35, 0.5, 0.65, 0.8):
        finite_difference = (outgoing_sigma(x + eps) - outgoing_sigma(x - eps)) / (2.0 * eps)
        assert outgoing_sigma_derivative(x) == pytest.approx(finite_difference, rel=2e-8, abs=2e-9)


def test_explicit_step_bound_is_a_conservative_global_witness() -> None:
    # The proof of 32 is analytic and documented in schedule_pressure.py.  This
    # dense sweep is only a regression smoke test against implementation drift.
    xs = np.linspace(1e-4, 1.0 - 1e-4, 10001)
    derivatives = np.array([outgoing_sigma_derivative(float(x)) for x in xs])
    assert np.all(derivatives >= -1e-15)
    assert float(np.max(derivatives)) < EXPLICIT_STEP_BOUND
    assert EXPLICIT_STEP_BOUND == 32.0
    assert EXPLICIT_STEP_BOUND >= 1.0
    assert EXPLICIT_FLATTEN_LENGTH == pytest.approx(330.0 * math.log(2.0) + 1.0)
    assert EXPLICIT_FLATTEN_LENGTH > 0.0


def test_schedule_shape_exponent_uses_constructive_flattening_schedule() -> None:
    endpoint = 3.25
    length = EXPLICIT_FLATTEN_LENGTH
    assert schedule_shape_exponent(endpoint - 7.0, endpoint) == 1.0
    assert schedule_shape_exponent(endpoint, endpoint) == 1.0
    assert schedule_shape_exponent(endpoint + 0.5 * length, endpoint) == pytest.approx(0.5, abs=2e-15)
    assert schedule_shape_exponent(endpoint + length, endpoint) == 0.0
    assert schedule_shape_exponent(endpoint + 2.0 * length, endpoint) == 0.0
    assert flatten_end(endpoint) == pytest.approx(endpoint + length)

    samples = [schedule_shape_exponent(endpoint + s * length, endpoint) for s in np.linspace(-0.2, 1.2, 101)]
    assert all(0.0 <= value <= 1.0 for value in samples)
    assert all(a >= b - 2e-15 for a, b in zip(samples, samples[1:]))


def test_schedule_pressure_scalar_inputs_fail_closed() -> None:
    with pytest.raises(ValueError):
        outgoing_sigma(float("nan"))
    with pytest.raises(ValueError):
        outgoing_sigma_derivative(float("inf"))
    with pytest.raises(ValueError):
        schedule_shape_exponent(0.0, 0.0, flatten_length=0.0)
    with pytest.raises(ValueError):
        flatten_end(0.0, flatten_length=-1.0)
