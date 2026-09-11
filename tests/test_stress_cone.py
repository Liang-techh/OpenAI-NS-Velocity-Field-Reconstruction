import math

import numpy as np
import pytest

from openai_ns_reconstruction.stress_cone import (
    PositiveSignedStressDecomposition,
    normal_magnitude,
    paper_ratio_margin,
    signed_covariance_matrix,
    squared_amplitude_coefficients,
    stress_target,
)


def test_explicit_squared_amplitudes_match_independent_linear_solve():
    a = 2.75
    b = 1.6
    scale_minus = 0.8
    scale_plus = 1.25
    m = 3.2
    t = 0.45

    matrix = signed_covariance_matrix(
        a=a, b=b, scale_minus=scale_minus, scale_plus=scale_plus
    )
    target = stress_target(m=m, t=t)
    explicit = squared_amplitude_coefficients(
        a=a,
        b=b,
        scale_minus=scale_minus,
        scale_plus=scale_plus,
        m=m,
        t=t,
    )

    # Independent numerical linear algebra cross-check.  Production does not use
    # np.linalg.solve, so this is not a re-evaluation of the same closed formula.
    independent = np.linalg.solve(matrix, target)
    assert explicit == pytest.approx(independent, rel=2e-14, abs=2e-14)


def test_positive_certificate_reconstructs_target_and_determinant():
    cert = PositiveSignedStressDecomposition(
        a=2.75,
        b=1.6,
        scale_minus=0.8,
        scale_plus=1.25,
        m=3.2,
        t=0.45,
    )

    assert cert.cone_margin > 0.0
    assert np.all(cert.squared_amplitudes > 0.0)
    assert np.all(cert.amplitudes > 0.0)
    assert np.square(cert.amplitudes) == pytest.approx(cert.squared_amplitudes)
    assert cert.reconstructed_target == pytest.approx(cert.target, rel=2e-14, abs=2e-14)

    independent_det = float(np.linalg.det(cert.matrix))
    assert cert.determinant == pytest.approx(independent_det, rel=2e-14, abs=2e-14)
    assert cert.determinant < 0.0


def test_paper_ratio_constructor_implies_reference_strict_cone():
    c_star = -0.5
    u_star = 2.0
    m = 3.0
    t = 0.5

    ratio_margin = paper_ratio_margin(c_star=c_star, u_star=u_star, m=m, t=t)
    assert ratio_margin == pytest.approx(
        u_star / math.sqrt(1.0 + u_star**2) - abs(c_star * t / m)
    )
    assert ratio_margin > 0.0

    cert = PositiveSignedStressDecomposition.from_paper_ratio(
        c_star=c_star,
        u_star=u_star,
        scale_minus=0.7,
        scale_plus=1.1,
        m=m,
        t=t,
    )
    assert cert.a == pytest.approx(-c_star * math.sqrt(1.0 + u_star**2))
    assert cert.b == pytest.approx(u_star)
    assert cert.cone_margin > 0.0
    assert cert.reconstructed_target == pytest.approx(cert.target, rel=2e-14, abs=2e-14)


def test_scaled_covariance_matches_epsilon_mask_squared_target():
    cert = PositiveSignedStressDecomposition(
        a=1.8,
        b=2.2,
        scale_minus=0.9,
        scale_plus=1.3,
        m=4.0,
        t=-0.6,
    )
    epsilon = 0.04
    mask = 0.7

    scaled = cert.scaled_reconstructed_target(epsilon=epsilon, mask=mask)
    expected = epsilon * mask**2 * cert.target
    assert scaled == pytest.approx(expected, rel=3e-14, abs=3e-14)

    # The sign of a finite mask disappears only through the covariance square.
    assert cert.scaled_reconstructed_target(epsilon=epsilon, mask=-mask) == pytest.approx(
        expected, rel=3e-14, abs=3e-14
    )


def test_strict_cone_boundary_and_outside_fail_closed():
    # |a*t| = b*m exactly: the reference theorem requires strict inequality.
    with pytest.raises(ValueError, match="strict cone"):
        PositiveSignedStressDecomposition(
            a=2.0,
            b=1.0,
            scale_minus=1.0,
            scale_plus=1.0,
            m=1.0,
            t=0.5,
        )

    with pytest.raises(ValueError, match="strict cone"):
        PositiveSignedStressDecomposition(
            a=2.0,
            b=1.0,
            scale_minus=1.0,
            scale_plus=1.0,
            m=1.0,
            t=0.6,
        )


def test_paper_ratio_boundary_fails_closed_instead_of_rounding_into_cone():
    c_star = -0.5
    u_star = 2.0
    m = 3.0
    boundary_t = m * u_star / (abs(c_star) * math.sqrt(1.0 + u_star**2))

    # Computing the boundary in binary64 may land a few ulps to either side, so
    # use the next representable value outside to make the rejected side strict.
    t_outside = math.nextafter(boundary_t, math.inf)
    assert paper_ratio_margin(c_star=c_star, u_star=u_star, m=m, t=t_outside) <= 0.0
    with pytest.raises(ValueError, match="paper reference ratio"):
        PositiveSignedStressDecomposition.from_paper_ratio(
            c_star=c_star,
            u_star=u_star,
            scale_minus=1.0,
            scale_plus=1.0,
            m=m,
            t=t_outside,
        )


def test_invalid_or_nonfinite_reference_inputs_are_rejected():
    with pytest.raises(ValueError, match="c_star<0"):
        normal_magnitude(0.1, 1.0)
    with pytest.raises(ValueError, match="scale_minus must be positive"):
        PositiveSignedStressDecomposition(
            a=1.0,
            b=1.0,
            scale_minus=0.0,
            scale_plus=1.0,
            m=2.0,
            t=0.0,
        )
    with pytest.raises(ValueError, match="finite"):
        squared_amplitude_coefficients(
            a=1.0,
            b=1.0,
            scale_minus=1.0,
            scale_plus=1.0,
            m=2.0,
            t=math.inf,
        )

    cert = PositiveSignedStressDecomposition(
        a=1.0,
        b=1.0,
        scale_minus=1.0,
        scale_plus=1.0,
        m=2.0,
        t=0.0,
    )
    with pytest.raises(ValueError, match="epsilon must be nonnegative"):
        cert.scaled_wave_amplitudes(epsilon=-1e-6, mask=1.0)
