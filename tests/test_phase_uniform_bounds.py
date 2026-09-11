import math

import pytest

from openai_ns_reconstruction.phase_estimates import PhaseScaleCertificate
from openai_ns_reconstruction.phase_uniform_bounds import (
    LocalBaseNumericEnvelope,
    UniformLocalBaseBridge,
)


def _envelope(**overrides):
    data = dict(
        M=4.0,
        epsilon=0.01,
        second_F0_bound=3.0,
        second_G0_bound=2.0,
        first_F0_bound=3.5,
        value_F_error_bound=3.0e-4,
        first_F_error_bound=2.0e-4,
        first_G_error_bound=4.0e-4,
        radial_F_bound=3.8,
        axial_F_bound=2.5,
        axial_G_bound=3.9,
    )
    data.update(overrides)
    return LocalBaseNumericEnvelope(**data)


def _scale(**overrides):
    data = dict(M=4.0, S=10.0, epsilon=0.01, k=100.0)
    data.update(overrides)
    return PhaseScaleCertificate(**data)


def test_uniform_bridge_matches_pinned_localbase_arithmetic():
    bridge = UniformLocalBaseBridge(local=_envelope(), scale=_scale(), diameter=8.0e-4)

    assert bridge.max_slow_box_diameter == pytest.approx(1.0e-3)
    assert bridge.derived_reference_error == pytest.approx(4.0 * (8.0e-4 + 1.0e-4))
    assert bridge.theorem_reference_error == pytest.approx(4.0 * (1.0e-3 + 1.0e-4))
    assert all(bridge.implication_checks().values())

    bounds = bridge.rounded_normal_local_bounds()
    assert bounds == pytest.approx(
        {
            "FR_abs": 4.0,
            "FZ_abs": 4.0,
            "GZ_abs": 4.0,
            "FR_minus_FR0_abs": 4.4e-3,
            "GR_minus_GR0_abs": 4.4e-3,
        }
    )


def test_envelope_fails_closed_when_localbase_threshold_is_exceeded():
    with pytest.raises(ValueError, match="second_F0_bound"):
        _envelope(second_F0_bound=4.000001)

    # M * epsilon^2 = 4e-4: the C1 actual/reference error may not exceed it.
    with pytest.raises(ValueError, match="first_G_error_bound"):
        _envelope(first_G_error_bound=4.000001e-4)


def test_bridge_fails_closed_outside_slow_box_diameter():
    with pytest.raises(ValueError, match="diameter"):
        UniformLocalBaseBridge(local=_envelope(), scale=_scale(), diameter=1.000001e-3)


def test_bridge_requires_exact_scale_parameter_alignment():
    with pytest.raises(ValueError, match="M must match"):
        UniformLocalBaseBridge(local=_envelope(), scale=_scale(M=5.0), diameter=5.0e-4)

    # This scale is still admissible: S=8, eps=.008 gives eps*S^2=.512.
    other = PhaseScaleCertificate(M=4.0, S=8.0, epsilon=0.008, k=64.0)
    with pytest.raises(ValueError, match="epsilon must match"):
        UniformLocalBaseBridge(local=_envelope(), scale=other, diameter=5.0e-4)


def test_independent_polynomial_derivative_family_obeys_derived_uniform_bound():
    """Independent analytic-family check; these polynomials are not paper fields."""
    local = _envelope()
    diameter = 8.0e-4
    bridge = UniformLocalBaseBridge(local=local, scale=_scale(), diameter=diameter)

    # Reference radial derivatives are affine with slopes bounded by the declared
    # second-derivative caps. Actual derivatives add independent C1 perturbations
    # no larger than M*eps^2. Compare every point directly with q0=0 rather than
    # using the bridge formula to manufacture the point values.
    delta_F = 2.0e-4
    delta_G = -3.0e-4
    slope_F = 3.0
    slope_G = -2.0
    q0_FR = 0.0
    q0_GR = 0.0

    for q in (-diameter, -0.4 * diameter, 0.0, 0.3 * diameter, diameter):
        actual_FR = slope_F * q + delta_F
        actual_GR = slope_G * q + delta_G
        assert abs(actual_FR - q0_FR) <= bridge.derived_reference_error
        assert abs(actual_GR - q0_GR) <= bridge.derived_reference_error

    assert math.isclose(local.first_G_error_bound, local.epsilon_squared_error_threshold)
