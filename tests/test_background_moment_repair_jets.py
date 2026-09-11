import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_jets import (
    solve_lemma52_eta_jets,
)


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


def _raw_derivative_jets(polynomial_coefficients):
    coefficients = np.asarray(polynomial_coefficients, dtype=float)
    factors = np.array([math.factorial(k) for k in range(coefficients.shape[0])])
    return coefficients * factors.reshape((-1,) + (1,) * (coefficients.ndim - 1))


def _trim_product(left, right, count):
    return np.polynomial.polynomial.polymul(left, right)[:count]


def test_eta_jet_solver_recovers_independently_constructed_polynomial_coefficients():
    repair = _repair()
    count = 4

    # Polynomial coefficients (not raw derivatives) for p(eta), alpha(eta),
    # beta(eta).  Base moments are generated in the *forward* Eq. (5.16)
    # direction using numpy polynomial multiplication, independently of the
    # quotient/Leibniz recurrences used by production.
    patch_coeff = np.array([1.8, -0.2, 0.05, 0.01])
    alpha_coeff = np.array(
        [[0.4, -0.2], [0.1, 0.3], [-0.04, 0.02], [0.01, -0.03]]
    )
    beta_coeff = np.array(
        [
            [0.15, -0.12, 0.08],
            [-0.05, 0.09, 0.04],
            [0.02, -0.03, 0.01],
            [-0.004, 0.006, -0.002],
        ]
    )

    u_first = alpha_coeff @ repair.u_matrix[0]
    u_weighted = alpha_coeff @ repair.u_matrix[1]
    e_first = beta_coeff @ repair.e_matrix[0]
    e_weighted = beta_coeff @ repair.e_matrix[1]
    e_last = beta_coeff @ repair.e_matrix[2]

    moment_coeff = np.zeros((count, 5))
    moment_coeff[:, 0] = -u_first
    moment_coeff[:, 1] = -e_first
    moment_coeff[:, 2] = -2.0 * _trim_product(patch_coeff, e_weighted, count)
    moment_coeff[:, 3] = -_trim_product(patch_coeff, u_weighted, count)
    moment_coeff[:, 4] = _trim_product(patch_coeff, e_last, count)

    out = solve_lemma52_eta_jets(
        repair,
        _raw_derivative_jets(moment_coeff),
        _raw_derivative_jets(patch_coeff),
    )

    assert out.max_order == 3
    assert np.allclose(
        out.alpha_derivatives,
        _raw_derivative_jets(alpha_coeff),
        rtol=3e-11,
        atol=3e-11,
    )
    assert np.allclose(
        out.beta_derivatives,
        _raw_derivative_jets(beta_coeff),
        rtol=3e-11,
        atol=3e-11,
    )
    assert np.max(np.abs(out.corrected_moment_derivatives)) < 2e-10


def test_order_zero_eta_jet_agrees_with_landed_pointwise_repair():
    repair = _repair()
    base = np.array([0.31, -0.27, 0.19, 0.41, -0.23])
    patch = 1.7
    pointwise = repair.solve(base, patch)
    jet = solve_lemma52_eta_jets(repair, base[None, :], np.array([patch]))

    assert jet.max_order == 0
    assert np.allclose(jet.alpha_derivatives[0], pointwise.alpha, rtol=2e-13, atol=2e-13)
    assert np.allclose(jet.beta_derivatives[0], pointwise.beta, rtol=2e-13, atol=2e-13)
    assert np.allclose(
        jet.corrected_moment_derivatives[0],
        pointwise.corrected_moments,
        rtol=0.0,
        atol=2e-13,
    )


def test_eta_jet_outputs_are_read_only():
    repair = _repair()
    out = solve_lemma52_eta_jets(
        repair,
        np.zeros((2, 5)),
        np.array([2.0, 0.1]),
    )
    with pytest.raises(ValueError):
        out.alpha_derivatives[0, 0] = 1.0
    with pytest.raises(ValueError):
        out.beta_derivatives[0, 0] = 1.0
    with pytest.raises(ValueError):
        out.corrected_moment_derivatives[0, 0] = 1.0


@pytest.mark.parametrize(
    "moments,patch,error",
    [
        (np.zeros(5), np.array([1.0]), ValueError),
        (np.zeros((2, 4)), np.array([1.0, 0.0]), ValueError),
        (np.zeros((2, 5)), np.array([1.0]), ValueError),
        (np.zeros((2, 5)), np.array([0.0, 1.0]), ValueError),
        (np.zeros((2, 5)), np.array([math.nan, 0.0]), ValueError),
    ],
)
def test_bad_eta_jet_data_fails_closed(moments, patch, error):
    with pytest.raises(error):
        solve_lemma52_eta_jets(_repair(), moments, patch)


def test_non_repair_object_fails_closed():
    with pytest.raises(TypeError):
        solve_lemma52_eta_jets(object(), np.zeros((1, 5)), np.ones(1))
