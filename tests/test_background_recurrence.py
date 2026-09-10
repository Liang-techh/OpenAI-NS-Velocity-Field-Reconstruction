import math

import pytest

from openai_ns_reconstruction.background_recurrence import (
    pressure_derivative_eq_5_5,
    pressure_row_terms_eq_5_5,
    solve_pressure_eq_5_5,
)


def test_order_one_pressure_row_has_the_exact_lemma_5_1_solution():
    """A polynomial n=1 instance checks Eq. (5.5), not a fitted surrogate.

    With C=2,
      phi_0 = 1 + eta + X,
      phi_1 = X (2-eta),
      Omega_0/X = 2(1+eta)X,
    Eq. (5.5) has an elementary primitive satisfying Pi_1(0,eta)=0.
    The source is specified independently of the implementation and the
    expected primitive below is derived by hand from the displayed equation.
    """

    C = 2.0
    phi = (
        lambda X, eta: 1.0 + eta + X,
        lambda X, eta: X * (2.0 - eta),
    )
    omega0_over_x = lambda X, eta: 2.0 * (1.0 + eta) * X

    for eta in (-0.8, -0.1, 0.4, 1.0):
        a = 0.5 * (2.0 - eta)
        for X in (0.0, 0.07, 0.4, 1.3):
            expected = 0.5 * (1.0 + eta) * (a - 1.0) * X**2 + a * X**3 / 3.0
            got = solve_pressure_eq_5_5(
                1, X, eta, phi, omega0_over_x, C=C, quadrature_points=8
            )
            assert math.isclose(got, expected, rel_tol=2e-14, abs_tol=2e-14)

    # Lemma 5.1 fixes every positive-order pressure trace at the axis.
    assert solve_pressure_eq_5_5(1, 0.0, 0.3, phi, omega0_over_x, C=C) == 0.0


def test_pressure_row_source_split_matches_the_full_convolution():
    """The post-(5.6) split must equal C^-2 sum_{i+j=n} phi_i phi_j."""

    phi = tuple(
        (lambda k: (lambda X, eta: (k + 1.0) + (k + 0.5) * X + eta))(k)
        for k in range(4)
    )
    omega_over_x = lambda X, eta: 0.7 - 0.2 * X + eta
    n, X, eta, C = 3, 0.6, -0.25, 2.3

    terms = pressure_row_terms_eq_5_5(n, X, eta, phi, omega_over_x, C=C)
    values = [f(X, eta) for f in phi]
    full = sum(values[i] * values[n - i] for i in range(n + 1)) / C**2
    full -= omega_over_x(X, eta) / 2.0

    assert math.isclose(terms.total, full, rel_tol=2e-15, abs_tol=2e-15)
    assert terms.lower_order_convolution != 0.0


def test_regular_omega_quotient_is_used_directly_at_the_axis():
    calls = []

    def quotient(X, eta):
        calls.append(X)
        return 3.0 + eta

    phi = (lambda X, eta: 2.0, lambda X, eta: X)
    derivative = pressure_derivative_eq_5_5(1, 0.0, 0.2, phi, quotient, C=4.0)
    assert math.isfinite(derivative)
    assert derivative == pytest.approx(-1.6)
    assert calls == [0.0]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n": 0},
        {"n": True},
        {"X": -1.0},
        {"X": float("nan")},
        {"eta": 1.1},
        {"C": 0.0},
        {"C": float("inf")},
    ],
)
def test_pressure_row_rejects_invalid_recurrence_data(kwargs):
    args = dict(
        n=1,
        X=0.2,
        eta=0.1,
        phi=(lambda X, eta: 1.0, lambda X, eta: X),
        omega_prev_over_x=lambda X, eta: 0.0,
        C=2.0,
    )
    args.update(kwargs)
    with pytest.raises(ValueError):
        pressure_derivative_eq_5_5(**args)


def test_pressure_row_requires_current_coefficient_and_regular_source():
    with pytest.raises(ValueError):
        pressure_derivative_eq_5_5(
            2,
            0.2,
            0.1,
            (lambda X, eta: 1.0, lambda X, eta: X),
            lambda X, eta: 0.0,
            C=2.0,
        )
    with pytest.raises(TypeError):
        pressure_derivative_eq_5_5(
            1,
            0.2,
            0.1,
            (lambda X, eta: 1.0, lambda X, eta: X),
            None,
            C=2.0,
        )
