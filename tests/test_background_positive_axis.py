import numpy as np
import pytest

from openai_ns_reconstruction.background_positive_axis import (
    PositiveAxisBaseJet,
    PositiveAxisSourceJet,
    RadialParameterJet,
    first_positive_order_term_from_jets,
    positive_axis_A0,
    positive_axis_A1,
    positive_axis_eq_5_7_fields,
    positive_axis_forcing,
)


def _data():
    base = PositiveAxisBaseJet(
        phi=RadialParameterJet(0.8, -0.3, 0.2, 0.11),
        axial=RadialParameterJet(-0.25, 0.4, -0.15, -0.07),
        beta=0.18,
    )
    source = PositiveAxisSourceJet(0.23, -0.31, 0.17, -0.09)
    return base, source


def test_matrix_rhs_matches_independent_expanded_system():
    h, order, C, xi, eta = 0.005, 1, 1.7, 0.41, 0.27
    lam = 2 * order * h
    X = xi * xi
    base, source = _data()
    phi = RadialParameterJet(0.14, -0.12, 0.05, 0.08)
    u = RadialParameterJet(-0.22, 0.09, -0.04, -0.06)
    k = RadialParameterJet(0.07, -0.03, 0.02, 0.05)
    p = RadialParameterJet(-0.13, 999.0, 0.01, 0.04)  # radial is replaced in the axial row
    q4, q5 = 0.37, -0.28

    w = np.array([
        phi.value, u.value, k.value, p.value,
        2 * xi * phi.radial, 2 * xi * u.radial,
    ])
    v = np.array([phi.parameter, u.parameter, k.parameter, p.parameter, q4, q5])
    got = (
        positive_axis_A0(h, order, C, xi, eta, base) @ w
        + positive_axis_A1(h, xi, eta, base) @ v
        + positive_axis_forcing(h, C, xi, eta, source)
    )

    a = 0.5 + h
    d = 0.5 - h
    edge = 1 - eta**2
    ell = 1 - 2 * h * eta**2
    angular_power = -a - 0.5
    axial_power = -a
    inv_c2 = 1 / C**2

    def time_value(power, jet):
        return (-power * jet.value + d * eta * jet.parameter + X * jet.radial) / ell

    def axial_value(power, jet):
        return (
            2 * eta * power * jet.value + edge * jet.parameter - 2 * eta * X * jet.radial
        ) / ell

    beta_value = (
        2 * eta * (a - lam) * u.value
        - 2 * eta * (d + lam) * k.value
        - edge * (u.parameter + k.parameter)
    ) / ell
    psrc = inv_c2 * source.pressure_product - source.omega_quotient / 2
    p_radial = 2 * inv_c2 * base.phi.value * phi.value + psrc
    p_sub = RadialParameterJet(p.value, p_radial, p.radial2, p.parameter)

    angular_rhs = (
        time_value(angular_power + lam, phi)
        + base.beta * (X * phi.radial + phi.value)
        + beta_value * (X * base.phi.radial + base.phi.value)
        + base.axial.value * axial_value(angular_power + lam, phi)
        + u.value * axial_value(angular_power, base.phi)
        + source.angular
    )
    axial_rhs = (
        time_value(axial_power + lam, u)
        + base.beta * X * u.radial
        + beta_value * X * base.axial.radial
        + base.axial.value * axial_value(axial_power + lam, u)
        + u.value * axial_value(axial_power, base.axial)
        + axial_value(-2 * a + lam, p_sub)
        + source.axial
    )
    expected = np.array([
        2 * xi * phi.radial,
        2 * xi * u.radial,
        -2 * xi * u.radial,
        2 * xi * p_radial,
        2 * angular_rhs,
        2 * axial_rhs,
    ])
    np.testing.assert_allclose(got, expected, rtol=5e-14, atol=5e-14)


def test_A1_has_the_pinned_sparse_shape_and_ignores_high_parameter_slots():
    base, _ = _data()
    A1 = positive_axis_A1(0.005, 0.4, -0.2, base)
    assert np.count_nonzero(A1[:4, :]) == 0
    assert np.count_nonzero(A1[:, 4:]) == 0
    v = np.array([0.1, -0.2, 0.3, -0.4, 900.0, -700.0])
    v2 = v.copy()
    v2[4:] = [-11.0, 13.0]
    np.testing.assert_array_equal(A1 @ v, A1 @ v2)


def test_first_positive_order_term_uses_exact_forcing_and_closed_form_G_integrals():
    h, C, eta, xi = 0.005, 1.5, 0.25, 0.6
    base, _ = _data()
    source = PositiveAxisSourceJet(
        angular=0.7,
        axial=-0.4,
        pressure_product=0.3,
        omega_quotient=-0.2,
    )
    got = first_positive_order_term_from_jets(
        1,
        xi,
        eta,
        h=h,
        C=C,
        base_provider=lambda _xi, _eta: base,
        source_provider=lambda _xi, _eta: source,
        quadrature_points=24,
    )
    psrc = source.pressure_product / C**2 - source.omega_quotient / 2
    ell = 1 - 2 * h * eta**2
    expected = np.array([
        0.0,
        0.0,
        0.0,
        psrc * xi**2,
        source.angular * xi / 2,
        source.axial * xi - eta * psrc * xi**3 / ell,
    ])
    np.testing.assert_allclose(got, expected, rtol=2e-13, atol=2e-14)


def test_provider_adapter_and_validation_fail_closed():
    base, source = _data()
    fields = positive_axis_eq_5_7_fields(
        0.005,
        2,
        2.0,
        lambda _xi, _eta: base,
        lambda _xi, _eta: source,
    )
    assert fields.A0(0.0, 0.0).shape == (6, 6)
    assert fields.A1(0.0, 0.0).shape == (6, 6)
    assert fields.forcing(0.0, 0.0).shape == (6,)
    with pytest.raises(ValueError):
        positive_axis_A0(0.005, 0, 1.0, 0.1, 0.0, base)
    with pytest.raises(ValueError):
        positive_axis_A1(0.005, 0.1, 1.1, base)
    with pytest.raises(ValueError):
        positive_axis_forcing(0.005, 0.0, 0.1, 0.0, source)
    bad = PositiveAxisSourceJet(0.0, 0.0, 0.0, 0.0)
    with pytest.raises(TypeError):
        positive_axis_eq_5_7_fields(
            0.005, 1, 1.0, lambda *_: base, lambda *_: (bad.angular,)
        ).forcing(0.1, 0.0)
