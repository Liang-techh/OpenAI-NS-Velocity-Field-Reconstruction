from fractions import Fraction
import pytest

from openai_ns_reconstruction.kokuno_profile_swirl_pressure import (
    FULL_RECONSTRUCTION,
    IMPORTED_PROFILE_EXISTENCE_PROVED_HERE,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    cartesian_swirl_energy_residual,
    cartesian_tangential_profile,
    pressure_gradient_from_swirl,
    profile_X_from_cartesian,
    regularized_swirl_energy,
    squared_norm2,
    swirl_pressure_residual,
)


def test_exact_swirl_pressure_identity():
    X = Fraction(13, 17)
    F = Fraction(-5, 7)
    pi_x = pressure_gradient_from_swirl(F)
    assert pi_x == Fraction(25, 49)
    assert regularized_swirl_energy(X, F) == Fraction(650, 833)
    assert swirl_pressure_residual(X, F, pi_x) == 0


def test_tiny_pressure_gradient_mutation_is_visible_exactly():
    X = Fraction(13, 17)
    F = Fraction(-5, 7)
    mutation = Fraction(1, 2**40)
    residual = swirl_pressure_residual(X, F, pressure_gradient_from_swirl(F) + mutation)
    assert residual == 2 * X * mutation
    assert residual != 0


def test_cartesian_axis_regular_swirl_energy_identity():
    y1, y2, F = Fraction(3, 5), Fraction(4, 5), Fraction(-7, 11)
    assert profile_X_from_cartesian(y1, y2) == Fraction(1, 2)
    vector = cartesian_tangential_profile(y1, y2, F)
    assert vector == (Fraction(28, 55), Fraction(-21, 55))
    assert squared_norm2(vector) == Fraction(49, 121)
    assert cartesian_swirl_energy_residual(y1, y2, F) == 0


def test_axis_tangential_factor_is_exactly_zero_for_finite_F():
    assert cartesian_tangential_profile(0, 0, Fraction(123, 17)) == (0, 0)
    assert cartesian_swirl_energy_residual(0, 0, Fraction(123, 17)) == 0


def test_rejects_inexact_or_nonphysical_inputs():
    with pytest.raises(TypeError):
        regularized_swirl_energy(0.5, Fraction(1, 2))
    with pytest.raises(TypeError):
        pressure_gradient_from_swirl(0.5)
    with pytest.raises(ValueError):
        regularized_swirl_energy(Fraction(-1, 10), Fraction(1, 2))


def test_target_leading_profile_swirl_pressure_interface():
    from openai_ns_reconstruction.profiles import LeadingProfile

    profile = LeadingProfile(
        E=lambda X, eta: (2 * X) ** 0.5 * (1.0 + eta / 8.0),
        F=lambda X, eta: 1.0 + eta / 8.0,
        U=lambda X, eta: 0.0,
        dU_deta=lambda X, eta: 0.0,
        Pi=lambda X, eta: X * (1.0 + eta / 8.0) ** 2,
        provenance="analytic interface regression only",
    )
    X, eta = 0.375, -0.5
    F = profile.smooth_swirl_factor(X, eta)
    assert F == pytest.approx(0.9375, rel=0.0, abs=0.0)
    assert profile.pressure_radial_derivative(X, eta) == pytest.approx(F * F, rel=0.0, abs=0.0)
    assert profile.paper_exact is False
    assert not PAPER_EXACT_VELOCITY_AVAILABLE
    assert not FULL_RECONSTRUCTION
    assert not IMPORTED_PROFILE_EXISTENCE_PROVED_HERE
