import math

import numpy as np

from openai_ns_reconstruction import LeadingProfile
from openai_ns_reconstruction.background import (
    BackgroundCoefficient,
    coefficient_radial_flux,
    coefficient_streamfunction,
    coefficient_vector_potential_cartesian,
)
from openai_ns_reconstruction.coordinates import similarity_coordinates


def polynomial_profile() -> LeadingProfile:
    # Deliberately simple coefficient data for kinematic identity tests only.
    # The tests do not claim this is one of the paper's recursively solved profiles.
    def E(X, eta):
        return 0.0

    def U(X, eta):
        return 1.0 + 0.4 * X + 0.3 * eta + 0.2 * eta * eta

    def dU_deta(X, eta):
        return 0.3 + 0.4 * eta

    return LeadingProfile(
        E=E,
        U=U,
        dU_deta=dU_deta,
        name="section5-kinematic-test-only",
        paper_exact=False,
    )


def test_eq_5_2_reduces_to_leading_flux_at_order_zero():
    h = 0.007
    profile = polynomial_profile()
    coefficient = BackgroundCoefficient(n=0, profile=profile)
    X, eta = 0.61, 0.27

    got = coefficient_radial_flux(X, eta, coefficient, h=h)
    expected = profile.V0(X, eta, h)

    assert math.isclose(got, expected, rel_tol=1e-13, abs_tol=1e-13)


def test_eq_5_27_streamfunction_recovers_axial_and_radial_coefficients():
    """Numerically differentiate S_n and check the two identities in Eq. (5.27)."""
    h = 0.006
    profile = polynomial_profile()
    coefficient = BackgroundCoefficient(n=3, profile=profile)
    r, z, t = 0.43, 0.11, 0.32
    point = similarity_coordinates(r, z, t, h)
    lam = coefficient.lambda_n(h)

    # Eq. (5.27): u_z,n = partial_s S_n, s = r^2/2.
    s_phys = 0.5 * r * r
    ds = 1.0e-6
    r_plus = math.sqrt(2.0 * (s_phys + ds))
    r_minus = math.sqrt(2.0 * (s_phys - ds))
    dS_ds = (
        coefficient_streamfunction(r_plus, z, t, coefficient, h=h)
        - coefficient_streamfunction(r_minus, z, t, coefficient, h=h)
    ) / (2.0 * ds)
    expected_uz = point.q ** (-point.A + lam) * profile.U(point.X, point.eta)
    assert math.isclose(dS_ds, expected_uz, rel_tol=2e-7, abs_tol=2e-8)

    # Eq. (5.27): r u_r,n = -partial_z S_n = q^lambda_n V_n.
    dz = 1.0e-6
    minus_dS_dz = -(
        coefficient_streamfunction(r, z + dz, t, coefficient, h=h)
        - coefficient_streamfunction(r, z - dz, t, coefficient, h=h)
    ) / (2.0 * dz)
    Vn = coefficient_radial_flux(point.X, point.eta, coefficient, h=h)
    expected_rur = point.q**lam * Vn
    assert math.isclose(minus_dS_dz, expected_rur, rel_tol=5e-7, abs_tol=5e-8)


def test_eq_5_27_cartesian_vector_potential_formula():
    h = 0.005
    profile = polynomial_profile()
    coefficient = BackgroundCoefficient(n=2, profile=profile)
    x, y, z, t = 0.31, -0.17, 0.08, 0.41
    r = math.hypot(x, y)
    point = similarity_coordinates(r, z, t, h)
    lam = coefficient.lambda_n(h)

    got = coefficient_vector_potential_cartesian(x, y, z, t, coefficient, h=h)

    # From Eq. (5.27) and r^2 = 2 q X:
    # A_n = (1/2) q^(-A+lambda_n) [F_n/X] (-y,x,0),
    # with F_n/X = A_X(U_n).
    AU = profile.radial_average_U(point.X, point.eta)
    factor = 0.5 * point.q ** (-point.A + lam) * AU
    expected = factor * np.array([-y, x, 0.0])

    assert np.allclose(got, expected, rtol=1e-12, atol=1e-12)
