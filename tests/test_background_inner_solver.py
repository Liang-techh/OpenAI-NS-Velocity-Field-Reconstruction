import numpy as np
import pytest

from openai_ns_reconstruction.background_inner_solver import (
    EQ_5_7_SINGULAR_WEIGHTS,
    first_picard_term_eq_5_7,
    picard_map_eq_5_7,
    singular_inverse_eq_5_7,
)


def test_eq_5_7_singular_inverse_matches_independent_polynomial_primitive():
    """Check the displayed G kernel against a hand-integrated vector source."""

    coefficients = np.array([1.5, -2.0, 0.75, 3.0, -1.25, 2.5])
    powers = np.array([0, 1, 2, 3, 1, 2], dtype=int)

    def source(s, eta):
        # eta-dependence makes sure the parameter is passed through unchanged.
        return coefficients * (1.0 + 0.2 * eta) * s**powers

    xi, eta = 0.73, -0.35
    got = singular_inverse_eq_5_7(xi, eta, source, quadrature_points=16)
    weights = np.asarray(EQ_5_7_SINGULAR_WEIGHTS)
    expected = (
        coefficients
        * (1.0 + 0.2 * eta)
        * xi ** (powers + 1)
        / (powers + weights + 1)
    )
    np.testing.assert_allclose(got, expected, rtol=2e-14, atol=2e-14)

    # Lemma 5.1 imposes zero positive-order axis data; the implementation does
    # not evaluate a numerical 0/0 kernel there.
    np.testing.assert_array_equal(singular_inverse_eq_5_7(0.0, eta, source), np.zeros(6))


def test_first_positive_order_picard_term_is_exact_k0_term_G_f1():
    """This is a real n=1 Picard-series step, not a fitted coefficient profile."""

    amplitudes = np.array([2.0, 1.0, -3.0, 0.5, 4.0, -2.0])

    def f1(s, eta):
        return amplitudes * (1.0 + eta) * s

    xi, eta = 0.4, 0.25
    got = first_picard_term_eq_5_7(1, xi, eta, f1, quadrature_points=12)
    weights = np.asarray(EQ_5_7_SINGULAR_WEIGHTS)
    # Integral_0^xi (s/xi)^c a(1+eta)s ds
    expected = amplitudes * (1.0 + eta) * xi**2 / (weights + 2.0)
    np.testing.assert_allclose(got, expected, rtol=2e-14, atol=2e-14)


def test_picard_map_includes_nonzero_K_step_with_independent_closed_form_oracle():
    """Exercise K=G(A0+A1*d_eta), not merely the lower-order source term."""

    amplitudes = np.array([1.0, -0.5, 2.0, 3.0, -1.0, 0.25])
    diagonal = np.array([0.4, -0.2, 0.7, 0.5, -0.3, 1.2])

    def previous(s, eta):
        return amplitudes * (1.0 + eta) * s

    def previous_eta_derivative(s, eta):
        return amplitudes * s

    def A0(s, eta):
        return np.diag(diagonal)

    # A sparse parameter-derivative coupling is included to exercise the A1 row.
    coupling = np.zeros((6, 6))
    coupling[4, 0] = 0.6
    coupling[5, 1] = -0.8

    def A1(s, eta):
        return coupling

    def zero_source(s, eta):
        return np.zeros(6)

    xi, eta = 0.55, -0.2
    got = picard_map_eq_5_7(
        1,
        xi,
        eta,
        previous,
        previous_eta_derivative,
        A0,
        A1,
        zero_source,
        quadrature_points=16,
    )

    rhs_amplitude = diagonal * amplitudes * (1.0 + eta) + coupling @ amplitudes
    weights = np.asarray(EQ_5_7_SINGULAR_WEIGHTS)
    expected = rhs_amplitude * xi**2 / (weights + 2.0)
    np.testing.assert_allclose(got, expected, rtol=2e-14, atol=2e-14)


@pytest.mark.parametrize("xi", [-1.0, float("nan"), float("inf")])
def test_singular_inverse_rejects_invalid_radial_coordinate(xi):
    with pytest.raises(ValueError):
        singular_inverse_eq_5_7(xi, 0.0, lambda s, eta: np.zeros(6))


def test_eq_5_7_interfaces_fail_closed_on_bad_shapes_and_orders():
    with pytest.raises(ValueError):
        singular_inverse_eq_5_7(0.2, 0.0, lambda s, eta: np.zeros(5))
    with pytest.raises(ValueError):
        first_picard_term_eq_5_7(0, 0.2, 0.0, lambda s, eta: np.zeros(6))
    with pytest.raises(ValueError):
        picard_map_eq_5_7(
            1,
            0.2,
            0.0,
            lambda s, eta: np.zeros(6),
            lambda s, eta: np.zeros(6),
            lambda s, eta: np.zeros((5, 5)),
            lambda s, eta: np.zeros((6, 6)),
            lambda s, eta: np.zeros(6),
        )
