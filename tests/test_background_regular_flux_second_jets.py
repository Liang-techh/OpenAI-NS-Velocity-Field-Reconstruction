import math

import pytest

from openai_ns_reconstruction.background_regular_flux_second_jets import (
    AxialThirdMixedJet,
    regular_flux_second_jet_eq_5_2,
)


def _P(X):
    X = float(X)
    return 1.0 + 0.3 * X + 0.2 * X**2 + 0.05 * X**3


def _P1(X):
    X = float(X)
    return 0.3 + 0.4 * X + 0.15 * X**2


def _P2(X):
    return 0.4 + 0.3 * float(X)


def _Q(eta):
    eta = float(eta)
    return 1.0 + 0.4 * eta + 0.15 * eta**2 + 0.07 * eta**3


def _Q1(eta):
    eta = float(eta)
    return 0.4 + 0.3 * eta + 0.21 * eta**2


def _Q2(eta):
    return 0.3 + 0.42 * float(eta)


def _Q3(_eta):
    return 0.42


def _u_jet(X, eta):
    return AxialThirdMixedJet(
        value=_P(X) * _Q(eta),
        radial=_P1(X) * _Q(eta),
        radial2=_P2(X) * _Q(eta),
        parameter=_P(X) * _Q1(eta),
        radial_parameter=_P1(X) * _Q1(eta),
        parameter2=_P(X) * _Q2(eta),
        radial2_parameter=_P2(X) * _Q1(eta),
        radial_parameter2=_P1(X) * _Q2(eta),
        parameter3=_P(X) * _Q3(eta),
    )


def _average_P(X):
    X = float(X)
    # Integral_0^1 P(sX) ds, evaluated independently in closed form.
    return 1.0 + 0.3 * X / 2.0 + 0.2 * X**2 / 3.0 + 0.05 * X**3 / 4.0


def _beta_value(X, eta, *, h, order):
    X = float(X)
    eta = float(eta)
    coefficient = 0.5 - h + 2.0 * order * h
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    average = _average_P(X)
    return (
        2.0 * eta * _P(X) * _Q(eta)
        - 2.0 * coefficient * eta * average * _Q(eta)
        - d * average * _Q1(eta)
    ) / ell


def _first_x(fn, X, eta, step=2.0e-5):
    return (fn(X + step, eta) - fn(X - step, eta)) / (2.0 * step)


def _second_x(fn, X, eta, step=2.0e-4):
    return (
        -fn(X + 2.0 * step, eta)
        + 16.0 * fn(X + step, eta)
        - 30.0 * fn(X, eta)
        + 16.0 * fn(X - step, eta)
        - fn(X - 2.0 * step, eta)
    ) / (12.0 * step * step)


def _first_eta(fn, X, eta, step=2.0e-5):
    return (fn(X, eta + step) - fn(X, eta - step)) / (2.0 * step)


def _second_eta(fn, X, eta, step=2.0e-4):
    return (
        -fn(X, eta + 2.0 * step)
        + 16.0 * fn(X, eta + step)
        - 30.0 * fn(X, eta)
        + 16.0 * fn(X, eta - step)
        - fn(X, eta - 2.0 * step)
    ) / (12.0 * step * step)


def _mixed(fn, X, eta, hx=8.0e-5, he=8.0e-5):
    return (
        fn(X + hx, eta + he)
        - fn(X + hx, eta - he)
        - fn(X - hx, eta + he)
        + fn(X - hx, eta - he)
    ) / (4.0 * hx * he)


def test_eq_5_2_regular_flux_second_jet_matches_independent_closed_form_oracle():
    h = 0.005
    order = 2
    X = 0.7
    eta = 0.23
    jet = regular_flux_second_jet_eq_5_2(h, order, X, eta, _u_jet, quadrature_points=24)
    fn = lambda x, e: _beta_value(x, e, h=h, order=order)

    # The oracle differentiates the independently integrated closed-form Eq. (5.2)
    # value; it does not reuse the production average-jet/product-rule formulas.
    assert math.isclose(jet.value, fn(X, eta), rel_tol=0.0, abs_tol=2e-13)
    assert math.isclose(jet.radial, _first_x(fn, X, eta), rel_tol=2e-7, abs_tol=2e-8)
    assert math.isclose(jet.radial2, _second_x(fn, X, eta), rel_tol=2e-6, abs_tol=8e-7)
    assert math.isclose(jet.parameter, _first_eta(fn, X, eta), rel_tol=2e-7, abs_tol=2e-8)
    assert math.isclose(
        jet.radial_parameter, _mixed(fn, X, eta), rel_tol=3e-6, abs_tol=8e-7
    )
    assert math.isclose(
        jet.parameter2, _second_eta(fn, X, eta), rel_tol=2e-6, abs_tol=8e-7
    )


def test_eq_5_2_regular_flux_second_jet_is_axis_regular_without_dividing_by_x():
    h = 0.005
    order = 1
    eta = -0.31
    jet = regular_flux_second_jet_eq_5_2(h, order, 0.0, eta, _u_jet)
    expected = _beta_value(0.0, eta, h=h, order=order)
    assert math.isfinite(jet.value)
    assert math.isclose(jet.value, expected, rel_tol=0.0, abs_tol=2e-13)
    assert all(
        math.isfinite(value)
        for value in (
            jet.radial,
            jet.radial2,
            jet.parameter,
            jet.radial_parameter,
            jet.parameter2,
        )
    )


def test_eq_5_2_regular_flux_second_jet_exposes_required_third_mixed_eta_data():
    def incomplete(_X, _eta):
        return object()

    with pytest.raises(TypeError):
        regular_flux_second_jet_eq_5_2(0.005, 1, 0.5, 0.1, incomplete)

    with pytest.raises(ValueError):
        AxialThirdMixedJet(
            value=0.0,
            radial=0.0,
            radial2=0.0,
            parameter=0.0,
            radial_parameter=0.0,
            parameter2=0.0,
            radial2_parameter=0.0,
            radial_parameter2=0.0,
            parameter3=float("nan"),
        )


def test_eq_5_2_regular_flux_second_jet_fails_closed_on_invalid_order_and_domain():
    with pytest.raises(ValueError):
        regular_flux_second_jet_eq_5_2(0.005, -1, 0.5, 0.1, _u_jet)
    with pytest.raises(ValueError):
        regular_flux_second_jet_eq_5_2(0.005, 1, -0.5, 0.1, _u_jet)
    with pytest.raises(ValueError):
        regular_flux_second_jet_eq_5_2(0.005, 1, 0.5, 1.1, _u_jet)
