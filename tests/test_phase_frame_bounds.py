import math

import pytest

from openai_ns_reconstruction.phase_estimates import phase_constant as phase_estimates_constant
from openai_ns_reconstruction.phase_frame_bounds import (
    FrameDampingEnvelope,
    base_phase_constant,
    coordinate_constant,
    damping_constant,
    damping_denominator,
    damping_error_bound,
    frame_coordinate_error_bound,
    normal_lower,
)


def test_base_phase_constants_match_pinned_definitions():
    M = 2.0
    u = 1.5

    assert damping_denominator(u) == pytest.approx((1.0 + u * u) ** 1.5)
    assert normal_lower(M, u) == pytest.approx(
        math.sqrt((1.0 / M) / (4.0 * (1.0 + u * u) ** 1.5))
    )
    assert base_phase_constant(M) == pytest.approx(8.0 * phase_estimates_constant(2.0 * M))


def test_general_frame_and_damping_bounds_match_theorem_arithmetic():
    M = 3.0
    A = 2.5
    b = 0.4
    delta = 0.01
    eta = 0.03

    G = M + 2.0 + 2.0 * A
    E = eta + 8.0 * (1.0 + A) * delta / b
    assert frame_coordinate_error_bound(M=M, A=A, b=b, delta=delta, eta=eta) == pytest.approx(
        16.0 * G * G * (1.0 + G) * E
    )
    assert damping_error_bound(M=M, A=A, delta=delta) == pytest.approx(
        4.0 * M * (2.0 * A + 5.0) * delta
    )


def test_large_band_specialization_recovers_named_lean_constants():
    M = 2.0
    u = 1.0
    S = 1.0e6
    B = 1.0
    viscosity = 2.5

    cert = FrameDampingEnvelope.from_large_band(
        M=M, S=S, u=u, B=B, viscosity=viscosity
    )

    assert cert.A == pytest.approx(3.0 * M)
    assert cert.b == pytest.approx(normal_lower(M, u))
    assert cert.delta == pytest.approx(base_phase_constant(M) / S)
    assert cert.eta == pytest.approx(16.0 * M * M / S)
    assert cert.frame_error == pytest.approx(coordinate_constant(M, u) / S)
    assert cert.damping_error == pytest.approx(damping_constant(M) / S)
    assert all(cert.specialized_constant_checks(u=u, S=S).values())


def test_envelope_fails_closed_on_scalar_theorem_hypotheses():
    M = 2.0
    u = 1.0
    b = normal_lower(M, u)

    with pytest.raises(ValueError, match="b<=B"):
        FrameDampingEnvelope(
            M=M, A=1.0, b=b, B=0.5 * b, delta=0.0, eta=0.0, viscosity=1.0
        )

    with pytest.raises(ValueError, match="B<=M"):
        FrameDampingEnvelope(
            M=M, A=1.0, b=b, B=M + 0.1, delta=0.0, eta=0.0, viscosity=1.0
        )

    with pytest.raises(ValueError, match="<=B/2"):
        FrameDampingEnvelope(
            M=M, A=1.0, b=b, B=1.0, delta=0.500001, eta=0.0, viscosity=1.0
        )

    with pytest.raises(ValueError, match=r"viscosity in \[0,4\]"):
        FrameDampingEnvelope(
            M=M, A=1.0, b=b, B=1.0, delta=0.01, eta=0.0, viscosity=4.0001
        )


def test_large_band_constructor_refuses_uncertified_smallness():
    # For M=2 the pinned phase constant is large; a small S cannot satisfy
    # delta<=B/2 with B<=M.  The constructor must reject instead of silently
    # treating a finite-band numerical example as a theorem-certified band.
    with pytest.raises(ValueError, match="<=B/2"):
        FrameDampingEnvelope.from_large_band(M=2.0, S=100.0, u=1.0, B=1.0, viscosity=2.0)
