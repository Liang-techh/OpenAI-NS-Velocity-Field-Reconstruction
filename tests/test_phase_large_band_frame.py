import math
from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_frame import LargeBandPhaseFrameCertificate
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


def _independent_family_constants(M: float, u: float) -> tuple[float, float, float]:
    """Rebuild the pinned family constants without calling production helpers."""

    frequency = M + M * (M + M**4) + (M + M**4) + M**2 + 4.0
    phase_at_twice_frequency = 8.0 * (2.0 * frequency) ** 3 + 2.0 * (2.0 * frequency) ** 4
    family_phase = 8.0 * phase_at_twice_frequency
    normal_lower = math.sqrt((1.0 / M) / (4.0 * (1.0 + u * u) ** 1.5))

    A = 3.0 * M
    G = M + 2.0 + 2.0 * A
    coordinate = 16.0 * G * G * (1.0 + G) * (
        16.0 * M * M + 8.0 * (1.0 + A) * family_phase / normal_lower
    )
    damping = 4.0 * M * (6.0 * M + 5.0) * family_phase
    return family_phase, coordinate, damping


def test_combined_gate_has_independent_first_frame_smallness_boundary_for_fixture():
    """The scale gate passes earlier; the family delta gate first passes at 94119."""

    h = Fraction(1, 200)
    M = 2.0
    u = 1.0
    B = 1.0
    family_phase, _, _ = _independent_family_constants(M, u)

    # Independent closed-form boundary for delta=D/ell^2 <= B/2.
    assert family_phase / (94118**2) > B / 2.0
    assert family_phase / (94119**2) < B / 2.0

    # Both bands are already far beyond the exact phase-scale threshold 23204.
    below_scale = LargeBandPhaseScaleCertificate(ell=94118, h=h, M=M)
    assert all(below_scale.scale_hypotheses().values())
    with pytest.raises(ValueError, match="<=B/2"):
        LargeBandPhaseFrameCertificate(
            scale=below_scale,
            u=u,
            B=B,
            viscosity=2.5,
        )

    at_scale = LargeBandPhaseScaleCertificate(ell=94119, h=h, M=M)
    cert = LargeBandPhaseFrameCertificate(
        scale=at_scale,
        u=u,
        B=B,
        viscosity=2.5,
    )
    assert all(cert.scalar_checks().values())


def test_combined_bounds_match_independent_closed_form_constants():
    ell = 100000
    M = 2.0
    u = 1.0
    B = 1.0
    S = float(ell**2)
    scale = LargeBandPhaseScaleCertificate(ell=ell, h=Fraction(1, 200), M=M)
    cert = LargeBandPhaseFrameCertificate(scale=scale, u=u, B=B, viscosity=2.5)

    family_phase, coordinate, damping = _independent_family_constants(M, u)
    local_phase_constant = 8.0 * M**3 + 2.0 * M**4

    assert cert.ell == ell
    assert cert.S_star == ell**2
    assert cert.local_phase_normal_bound == pytest.approx(4.0 * local_phase_constant / S)
    assert cert.family_normal_delta == pytest.approx(family_phase / S)
    assert cert.frame_error_bound == pytest.approx(coordinate / S)
    assert cert.damping_error_bound == pytest.approx(damping / S)


def test_bridge_keeps_local_and_family_phase_constants_distinct():
    scale = LargeBandPhaseScaleCertificate(ell=100000, h=0.005, M=2.0)
    cert = LargeBandPhaseFrameCertificate(scale=scale, u=1.0, B=1.0, viscosity=1.0)

    assert cert.family_normal_delta > cert.local_phase_normal_bound
    assert cert.status == "formal-structure"
    assert cert.actual_base_fields_verified is False
    assert cert.local_base_vector_hypotheses_verified is False
    assert cert.uniform_eq_7_9_to_7_11_verified is False
    assert cert.paper_exact_velocity_available is False


def test_bridge_fails_closed_on_wrong_scale_type_and_frame_scalar_hypotheses():
    with pytest.raises(TypeError, match="LargeBandPhaseScaleCertificate"):
        LargeBandPhaseFrameCertificate(scale=object(), u=1.0, B=1.0, viscosity=1.0)

    scale = LargeBandPhaseScaleCertificate(ell=100000, h=Fraction(1, 200), M=2.0)

    with pytest.raises(ValueError, match="b<=B"):
        LargeBandPhaseFrameCertificate(scale=scale, u=1.0, B=0.1, viscosity=1.0)

    with pytest.raises(ValueError, match="B<=M"):
        LargeBandPhaseFrameCertificate(scale=scale, u=1.0, B=2.1, viscosity=1.0)

    with pytest.raises(ValueError, match=r"viscosity in \[0,4\]"):
        LargeBandPhaseFrameCertificate(scale=scale, u=1.0, B=1.0, viscosity=4.1)

    for bad in (float("nan"), float("inf"), True):
        with pytest.raises(ValueError):
            LargeBandPhaseFrameCertificate(scale=scale, u=bad, B=1.0, viscosity=1.0)
