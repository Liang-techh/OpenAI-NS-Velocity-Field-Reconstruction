import math
import numpy as np
import pytest

from openai_ns_reconstruction.phase_box_certificate import (
    RepresentativePhaseData,
    RoundedNormalPointCertificate,
    lean_nonzero_round,
    lean_rounded_frequency,
)
from openai_ns_reconstruction.phase_estimates import PhaseScaleCertificate


def test_pinned_lean_rounding_is_floor_with_zero_replaced_by_one():
    assert lean_nonzero_round(1.9) == 1
    assert lean_nonzero_round(-1.2) == -2
    assert lean_nonzero_round(0.2) == 1
    assert lean_nonzero_round(0.0) == 1

    # The pinned Lean quantitative theorem is not using nearest-integer rounding.
    assert lean_rounded_frequency(10.0, 0.19) == pytest.approx(0.1)
    assert abs(lean_rounded_frequency(10.0, 0.19) - 0.19) <= 0.1


def test_representative_frequency_equalities_hold_by_construction():
    rep = RepresentativePhaseData(
        R0=1.0,
        FR0=1.0,
        GR0=0.0,
        B=1.0,
        sigma=1,
        u=0.25,
        L=100.0,
    )
    residuals = rep.algebraic_residuals()
    assert abs(residuals["K_norm_minus_1"]) < 1e-15
    assert abs(residuals["K_dot_g"]) < 1e-15
    assert abs(residuals["frequency_theta"]) < 1e-15
    assert abs(residuals["frequency_z"]) < 1e-15

    expected = rep.B * (
        rep.K - rep.sigma * rep.u / (rep.L * rep.g_norm**2) * rep.g
    )
    assert np.allclose(rep.representative_frequency, expected, rtol=0.0, atol=1e-15)
    assert np.allclose(
        np.array([rep.target / rep.R0, rep.pz]),
        rep.representative_frequency,
        rtol=0.0,
        atol=1e-15,
    )


def test_point_certificate_instantiates_rounded_normal_estimate_hypotheses():
    rep = RepresentativePhaseData(
        R0=1.0,
        FR0=1.0,
        GR0=0.0,
        B=1.0,
        sigma=1,
        u=0.0,
        L=100.0,
    )
    scale = PhaseScaleCertificate(M=1.0, S=100.0, epsilon=1e-6, k=1e6)
    cert = RoundedNormalPointCertificate(
        representative=rep,
        scale=scale,
        R=1.0,
        v=0.0,
        FR=1.0,
        GR=0.0,
        FZ=0.0,
        GZ=0.0,
    )

    assert all(cert.hypotheses().values())
    assert cert.p == pytest.approx(1e-6)  # punctured lattice at target=0
    assert cert.normal_error == pytest.approx(1e-6)
    assert cert.normal_velocity_norm == pytest.approx(1e-6)
    assert cert.numeric_crosscheck() == {
        "normal_error_le_bound": True,
        "normal_velocity_le_bound": True,
    }
    assert cert.geometry_smallness()
    geometry = cert.geometry_crosscheck()
    assert geometry["normal_scale_ge_B_over_2"]
    assert geometry["normal_norm_ge_B_over_2"]
    assert geometry["radial_slope_error"] <= geometry["radial_slope_bound"]
    assert geometry["direction_error"] <= geometry["direction_bound"]


def test_point_certificate_fails_closed_on_missing_local_base_bound():
    rep = RepresentativePhaseData(
        R0=1.0,
        FR0=1.0,
        GR0=0.0,
        B=1.0,
        sigma=1,
        u=0.0,
        L=100.0,
    )
    scale = PhaseScaleCertificate(M=1.0, S=100.0, epsilon=1e-6, k=1e6)
    with pytest.raises(ValueError, match="FR_reference_error"):
        RoundedNormalPointCertificate(
            representative=rep,
            scale=scale,
            R=1.0,
            v=0.0,
            FR=1.01,
            GR=0.0,
            FZ=0.0,
            GZ=0.0,
        )


def test_representative_data_rejects_degenerate_shear_or_nonunit_signs():
    with pytest.raises(ValueError, match="shear g must be nonzero"):
        RepresentativePhaseData(
            R0=1.0, FR0=0.0, GR0=0.0, B=1.0, sigma=1, u=0.0, L=10.0
        )
    with pytest.raises(ValueError, match="sigma"):
        RepresentativePhaseData(
            R0=1.0, FR0=1.0, GR0=0.0, B=1.0, sigma=0, u=0.0, L=10.0
        )
    with pytest.raises(ValueError, match="orientation"):
        RepresentativePhaseData(
            R0=1.0,
            FR0=1.0,
            GR0=0.0,
            B=1.0,
            sigma=1,
            u=0.0,
            L=10.0,
            orientation=0,
        )


def test_negative_numeric_slack_is_rejected():
    rep = RepresentativePhaseData(
        R0=1.0, FR0=1.0, GR0=0.0, B=1.0, sigma=1, u=0.0, L=100.0
    )
    scale = PhaseScaleCertificate(M=1.0, S=100.0, epsilon=1e-6, k=1e6)
    cert = RoundedNormalPointCertificate(
        representative=rep,
        scale=scale,
        R=1.0,
        v=0.0,
        FR=1.0,
        GR=0.0,
        FZ=0.0,
        GZ=0.0,
    )
    with pytest.raises(ValueError):
        cert.numeric_crosscheck(slack=-1.0)
