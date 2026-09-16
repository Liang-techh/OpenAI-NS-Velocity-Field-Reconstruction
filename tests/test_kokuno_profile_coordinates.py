from fractions import Fraction as Q

import pytest

from openai_ns_reconstruction.kokuno_profile_coordinates import (
    FULL_RECONSTRUCTION,
    IMPORTED_PROFILE_EXISTENCE_PROVED_HERE,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    AxisymmetricFieldJet,
    ProfileJet,
    SelfSimilarPoint,
    incompressibility_profile_residual,
    omega_theta_over_r,
    scaled_profile_axial_derivative,
    scaled_profile_time_derivative,
)


def test_exact_coordinate_jacobian_and_inverse_factors() -> None:
    point = SelfSimilarPoint(h=Q(1, 100), q=Q(9, 4), eta=Q(2, 5), X=Q(7, 3))

    assert point.A == Q(51, 100)
    assert point.D == Q(49, 100)
    assert point.d == Q(21, 25)
    assert point.L == Q(623, 625)
    assert point.scaled_tau_z_jacobian == Q(623, 625)

    factors = point.inverse_derivative_factors()
    assert factors.q_t_factor == Q(-625, 623)
    assert factors.eta_t_factor == Q(245, 1246)
    assert factors.X_t_factor == Q(4375, 1869)
    assert factors.q_z_factor == Q(500, 623)
    assert factors.eta_z_factor == Q(525, 623)
    assert factors.X_z_factor == Q(-3500, 1869)


def test_exact_chain_rule_factors() -> None:
    point = SelfSimilarPoint(h=Q(1, 20), q=3, eta=Q(-1, 3), X=Q(5, 2))
    jet = ProfileJet(value=Q(7, 5), d_eta=Q(-11, 7), d_X=Q(13, 9))
    b = Q(-3, 4)

    expected_t = (
        -b * jet.value + point.D * point.eta * jet.d_eta + point.X * jet.d_X
    ) / point.L
    expected_z = (
        2 * b * point.eta * jet.value
        + point.d * jet.d_eta
        - 2 * point.eta * point.X * jet.d_X
    ) / point.L
    assert scaled_profile_time_derivative(point, b=b, jet=jet) == expected_t
    assert scaled_profile_axial_derivative(point, b=b, jet=jet) == expected_z


def test_incompressibility_is_exact_and_tiny_defect_does_not_pass() -> None:
    point = SelfSimilarPoint(h=Q(1, 200), q=2, eta=Q(3, 8), X=Q(4, 3))
    U = ProfileJet(value=Q(5, 9), d_eta=Q(-2, 7), d_X=Q(11, 13))
    certified_v0_x = (
        2 * point.A * point.eta * U.value
        - point.d * U.d_eta
        + 2 * point.eta * point.X * U.d_X
    ) / point.L

    assert incompressibility_profile_residual(point, V0_X=certified_v0_x, U=U) == 0
    assert incompressibility_profile_residual(
        point, V0_X=certified_v0_x + Q(1, 2**40), U=U
    ) == Q(1, 2**40)


def test_omega_theta_over_r_component_identity() -> None:
    # Polynomial witness: u_r=a*r*z and u_z=b*s, s=r^2/2. Then
    # d_z(r*u_r)=2*a*s and d_s(u_z)=b, hence omega_theta/r=a-b.
    a = Q(7, 11)
    b = Q(-5, 13)
    s = Q(9, 10)
    jet = AxisymmetricFieldJet(s=s, d_z_ru_r=2 * a * s, d_s_u_z=b)

    assert omega_theta_over_r(jet) == a - b


def test_fail_closed_domains_and_truth_boundary() -> None:
    with pytest.raises(ValueError, match=r"0 < h < 1/2"):
        SelfSimilarPoint(h=Q(1, 2), q=1, eta=0, X=0)
    with pytest.raises(ValueError, match=r"\|eta\| < 1"):
        SelfSimilarPoint(h=Q(1, 100), q=1, eta=1, X=0)
    with pytest.raises(TypeError, match="must be an int or Fraction"):
        SelfSimilarPoint(h=0.01, q=1, eta=0, X=0)
    with pytest.raises(ValueError, match=r"s=r\^2/2 > 0"):
        AxisymmetricFieldJet(s=0, d_z_ru_r=0, d_s_u_z=0)

    assert IMPORTED_PROFILE_EXISTENCE_PROVED_HERE is False
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
