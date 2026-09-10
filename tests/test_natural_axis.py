import math

import pytest

from openai_ns_reconstruction.natural_axis import (
    D,
    H,
    NaturalAxisParameters,
    NaturalProfileAssembly,
    W,
    axis_U,
    real_amplitude,
    real_phase,
)


def test_axis_polynomials_match_formalized_identities():
    h = 1.0e-4
    j = 2.0e-4
    eta = 0.37

    # NaturalAxisData.neg_W_formula.
    rhs = 3.0 - 8.0 * h * eta**2 + 2.0 * D(h) * j * eta
    assert -W(h, j, eta) == pytest.approx(rhs, rel=0.0, abs=2.0e-15)

    # NaturalAxisData.H_at_left / H_at_right.
    assert H(h, j, -j / 4.0) == pytest.approx(-(D(h) * j) / 4.0, abs=2.0e-20)
    right = (j / 5.0) * (1.0 - D(h) - j**2 / 25.0)
    assert H(h, j, -j / 5.0) == pytest.approx(right, abs=2.0e-20)


def test_real_phase_and_amplitude_at_axis_origin():
    h = 1.0e-4
    j = 2.0e-4
    sigma = 0.03
    Lambda = 7.0
    C = 2.5

    assert real_phase(h, j, sigma, 0.0) == 0.0
    assert real_amplitude(h, j, sigma, Lambda, C, 0.0) == pytest.approx(1.0 / C)


def test_natural_profile_assembly_matches_lean_rescaling_formulas():
    p = NaturalAxisParameters(
        h=1.0e-4,
        j=2.0e-4,
        sigma=0.03,
        Lambda=2.0,
        C=3.0,
        phase_samples=101,
    )

    phi = lambda Y, eta: Y + 2.0 * eta
    u = lambda Y, eta: Y * eta
    du_deta = lambda Y, eta: Y
    average = lambda Y, eta: Y**2
    pressure = lambda Y, eta: 3.0 * Y - eta
    axis_pressure = lambda eta: -1.0 - eta**2

    a = NaturalProfileAssembly(p, phi, u, du_deta, average, pressure, axis_pressure)
    X = 0.4
    eta = 0.2
    Y = p.Lambda * X

    assert a.rescale_point(X, eta) == pytest.approx((Y, eta))
    assert a.E(X, eta) == pytest.approx(p.amplitude(eta) * phi(Y, eta))
    assert a.U(X, eta) == pytest.approx(axis_U(p.j, eta) + u(Y, eta) / p.Lambda)
    assert a.dU_deta(X, eta) == pytest.approx(4.0 + du_deta(Y, eta) / p.Lambda)
    assert a.radial_average(X, eta) == pytest.approx(
        axis_U(p.j, eta) + average(Y, eta) / p.Lambda
    )
    assert a.Pi(X, eta) == pytest.approx(axis_pressure(eta) + pressure(Y, eta) / p.Lambda)

    leading = a.to_leading_profile()
    assert leading.paper_exact is False
    assert leading.E(X, eta) == pytest.approx(a.E(X, eta))
    assert leading.U(X, eta) == pytest.approx(a.U(X, eta))


def test_parameter_guard_matches_small_parameter_range():
    with pytest.raises(ValueError):
        NaturalAxisParameters(0.0, 1.0e-4, 0.03, 2.0, 3.0)
    with pytest.raises(ValueError):
        NaturalAxisParameters(1.0e-4, 2.0e-3, 0.03, 2.0, 3.0)

    p = NaturalAxisParameters(1.0e-3, 1.0e-3, 0.03, 2.0, 3.0)
    assert math.isfinite(p.amplitude(0.1))
