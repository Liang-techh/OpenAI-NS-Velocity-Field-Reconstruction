import math

import pytest
from scipy import special

from openai_ns_reconstruction.axis_analytic_input_bounds import (
    constructive_analytic_input_norm_certificate,
    factorial_resolvent_majorant_upper,
)


def test_factorial_resolvent_majorant_zero_is_exact() -> None:
    cert = factorial_resolvent_majorant_upper(0.0)
    assert cert.partial_term_count == 1
    assert cert.tail_ratio_upper == 0.0
    assert cert.series_upper == 1.0


def test_factorial_resolvent_majorant_crosschecks_bessel_identity() -> None:
    # Independent closed form:
    # sum K^k/(k!(k+1)!) = I_1(2 sqrt(K)) / sqrt(K).
    K = 30_720.0
    cert = factorial_resolvent_majorant_upper(K)
    exact_numeric = special.iv(1, 2.0 * math.sqrt(K)) / math.sqrt(K)

    assert cert.tail_ratio_upper <= 0.5
    assert cert.partial_term_count > 1
    assert math.isfinite(cert.series_upper)
    assert cert.series_upper >= exact_numeric * (1.0 - 1.0e-12)
    assert cert.series_upper / exact_numeric < 1.0 + 1.0e-10


def test_canonical_half_radius_collapses_norm_inputs() -> None:
    cert = constructive_analytic_input_norm_certificate(
        neighborhood_radius=0.2,
        complex_field_sup_upper=1.0,
    )

    assert cert.epsilon == 0.1
    assert cert.radius_loss_half == 12.0
    assert cert.amplitude_norm_upper == 12.0
    # The implementation rounds positive products upward, so the common family
    # norm bound is allowed to be one ulp above the exact 12*B value.
    assert cert.coefficient_family_norm_upper >= 12.0
    assert cert.chi_norm_upper == cert.coefficient_family_norm_upper
    assert cert.resolvent_majorant.majorant_parameter_upper >= 30_720.0

    data = cert.axis_data_norm_bounds(h=0.001)
    common = cert.coefficient_family_norm_upper
    assert data.A == pytest.approx(0.501)
    assert data.D == pytest.approx(0.499)
    assert data.h == 0.001
    assert data.one == common
    assert data.eta == common
    assert data.d == common
    assert data.inverse_l == common
    assert data.u_star == common
    assert data.u_star_eta == common
    assert data.w_star == common
    assert data.h_star == common
    assert data.normalized_gradient == common
    assert data.z_star == common


def test_canonical_certificate_instantiates_existing_operator_ledger() -> None:
    cert = constructive_analytic_input_norm_certificate(
        neighborhood_radius=2.0,
        complex_field_sup_upper=1.0,
    )
    operators = cert.operator_norm_bounds()

    assert operators.product == 64.0
    assert operators.average == 1.0
    assert operators.primitive == 80.0
    assert operators.parameter_primitive == 80.0
    assert operators.param1 == 5120.0
    assert operators.resolvent == cert.resolvent_majorant.series_upper


def test_canonical_certificate_fails_closed_on_missing_analytic_witnesses() -> None:
    with pytest.raises(ValueError, match="neighborhood_radius"):
        constructive_analytic_input_norm_certificate(
            neighborhood_radius=0.0,
            complex_field_sup_upper=1.0,
        )
    with pytest.raises(ValueError, match="complex_field_sup_upper"):
        constructive_analytic_input_norm_certificate(
            neighborhood_radius=1.0,
            complex_field_sup_upper=float("nan"),
        )
    with pytest.raises(ValueError, match="majorant_parameter_upper"):
        factorial_resolvent_majorant_upper(-1.0)
