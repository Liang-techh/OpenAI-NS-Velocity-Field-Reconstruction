import math

import pytest

from openai_ns_reconstruction.background_recurrence import (
    RegularFluxJet,
    omega_over_x_eq_5_6,
    omega_row_terms_eq_5_6,
    regular_T_on_Xv_over_X,
    regular_Z_on_Xv_over_X,
    regular_Z2_on_Xv_over_X,
)


def _jet0(X, eta):
    # v0 = 1 + .3 X - .2 eta + .4 X eta + .1 X^2 + .05 eta^2
    return RegularFluxJet(
        1.0 + 0.3 * X - 0.2 * eta + 0.4 * X * eta + 0.1 * X**2 + 0.05 * eta**2,
        0.3 + 0.4 * eta + 0.2 * X,
        -0.2 + 0.4 * X + 0.1 * eta,
        0.2,
        0.4,
        0.1,
    )


def _jet1(X, eta):
    # v1 = .6 - .1 X + .25 eta - .15 X eta + .07 X^2 + .03 eta^2
    return RegularFluxJet(
        0.6 - 0.1 * X + 0.25 * eta - 0.15 * X * eta + 0.07 * X**2 + 0.03 * eta**2,
        -0.1 - 0.15 * eta + 0.14 * X,
        0.25 - 0.15 * X + 0.06 * eta,
        0.14,
        -0.15,
        0.06,
    )


def _U0(X, eta):
    return 1.2 + 0.2 * X - 0.1 * eta


def _U1(X, eta):
    return -0.4 + 0.3 * X + 0.05 * eta


def _raw_V(jet_fn, X, eta):
    return X * jet_fn(X, eta).value


def _raw_V_X(jet_fn, X, eta):
    j = jet_fn(X, eta)
    return j.value + X * j.dX


def _raw_V_XX(jet_fn, X, eta):
    j = jet_fn(X, eta)
    return 2.0 * j.dX + X * j.dXX


def _raw_V_eta(jet_fn, X, eta):
    return X * jet_fn(X, eta).dEta


def _raw_T(jet_fn, X, eta, b, h):
    D = 0.5 - h
    L = 1.0 - 2.0 * h * eta**2
    V = _raw_V(jet_fn, X, eta)
    return (
        -b * V
        + D * eta * _raw_V_eta(jet_fn, X, eta)
        + X * _raw_V_X(jet_fn, X, eta)
    ) / L


def _raw_Z_from_values(value, dX, dEta, X, eta, b, h):
    d = 1.0 - eta**2
    L = 1.0 - 2.0 * h * eta**2
    return (2.0 * b * eta * value + d * dEta - 2.0 * eta * X * dX) / L


def _raw_Z(jet_fn, X, eta, b, h):
    return _raw_Z_from_values(
        _raw_V(jet_fn, X, eta),
        _raw_V_X(jet_fn, X, eta),
        _raw_V_eta(jet_fn, X, eta),
        X,
        eta,
        b,
        h,
    )


def _raw_Z2_fd(jet_fn, X, eta, b, h, eps=2e-6):
    """Independent finite-difference check of Z_{b-D}(Z_b V)."""

    D = 0.5 - h

    def first(x, e):
        return _raw_Z(jet_fn, x, e, b, h)

    dX = (first(X + eps, eta) - first(X - eps, eta)) / (2.0 * eps)
    dEta = (first(X, eta + eps) - first(X, eta - eps)) / (2.0 * eps)
    return _raw_Z_from_values(
        first(X, eta), dX, dEta, X, eta, b - D, h
    )


def _raw_omega_over_x(k, X, eta, h):
    """Direct Eq. (5.6) away from the axis, before using V=Xv regularity."""

    jets = (_jet0, _jet1)
    U = (_U0, _U1)
    lam = lambda n: 2.0 * n * h

    omega = _raw_T(jets[k], X, eta, lam(k), h)
    for i in range(k + 1):
        j = k - i
        V_i = _raw_V(jets[i], X, eta)
        V_j = _raw_V(jets[j], X, eta)
        omega += V_i * (_raw_V_X(jets[j], X, eta) - V_j / (2.0 * X))
        omega += U[i](X, eta) * _raw_Z(jets[j], X, eta, lam(j), h)
    omega -= 2.0 * X * _raw_V_XX(jets[k], X, eta)
    if k >= 1:
        omega -= _raw_Z2_fd(jets[k - 1], X, eta, lam(k - 1), h)
    return omega / X


def test_regular_T_and_Z_match_the_raw_eq_4_2_operators():
    X, eta, h = 0.43, -0.27, 0.008
    jet = _jet1(X, eta)
    b = 2.0 * h

    assert regular_T_on_Xv_over_X(jet, X, eta, b=b, h=h) == pytest.approx(
        _raw_T(_jet1, X, eta, b, h) / X, rel=2e-14, abs=2e-14
    )
    assert regular_Z_on_Xv_over_X(jet, X, eta, b=b, h=h) == pytest.approx(
        _raw_Z(_jet1, X, eta, b, h) / X, rel=2e-14, abs=2e-14
    )


def test_regular_shifted_Z2_matches_independent_raw_finite_difference():
    X, eta, h = 0.61, 0.23, 0.006
    jet = _jet0(X, eta)
    got = regular_Z2_on_Xv_over_X(jet, X, eta, b=0.0, h=h)
    expected = _raw_Z2_fd(_jet0, X, eta, 0.0, h) / X
    assert got == pytest.approx(expected, rel=3e-8, abs=3e-8)


@pytest.mark.parametrize("k", [0, 1])
def test_regular_omega_quotient_matches_direct_eq_5_6_away_from_axis(k):
    h, eta = 0.007, 0.31
    for X in (0.19, 0.57, 1.1):
        jets = (_jet0(X, eta), _jet1(X, eta))
        U = (_U0(X, eta), _U1(X, eta))
        got = omega_over_x_eq_5_6(k, X, eta, jets, U, h=h)
        expected = _raw_omega_over_x(k, X, eta, h)
        # k=1 contains an intentionally independent finite-difference Z^[2] check.
        tol = 5e-8 if k else 2e-14
        assert got == pytest.approx(expected, rel=tol, abs=tol)


def test_omega_quotient_is_regular_at_axis_and_uses_negative_index_zero_rule():
    h, eta = 0.005, -0.4
    jets0 = (_jet0(0.0, eta),)
    terms0 = omega_row_terms_eq_5_6(0, 0.0, eta, jets0, (_U0(0.0, eta),), h=h)
    assert math.isfinite(terms0.total)
    assert terms0.shifted_axial_viscosity == 0.0

    jets1_axis = (_jet0(0.0, eta), _jet1(0.0, eta))
    U_axis = (_U0(0.0, eta), _U1(0.0, eta))
    axis = omega_over_x_eq_5_6(1, 0.0, eta, jets1_axis, U_axis, h=h)
    assert math.isfinite(axis)

    X = 1e-7
    near = omega_over_x_eq_5_6(
        1,
        X,
        eta,
        (_jet0(X, eta), _jet1(X, eta)),
        (_U0(X, eta), _U1(X, eta)),
        h=h,
    )
    assert near == pytest.approx(axis, rel=2e-6, abs=2e-6)


def test_eq_5_6_rejects_incomplete_or_nonfinite_jets():
    good = RegularFluxJet(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    bad = RegularFluxJet(float("nan"), 0.0, 0.0, 0.0, 0.0, 0.0)

    with pytest.raises(ValueError):
        omega_over_x_eq_5_6(-1, 0.2, 0.0, (good,), (1.0,), h=0.005)
    with pytest.raises(ValueError):
        omega_over_x_eq_5_6(1, 0.2, 0.0, (good,), (1.0, 2.0), h=0.005)
    with pytest.raises(ValueError):
        omega_over_x_eq_5_6(0, 0.2, 0.0, (good,), (), h=0.005)
    with pytest.raises(ValueError):
        omega_over_x_eq_5_6(0, 0.2, 0.0, (bad,), (1.0,), h=0.005)
