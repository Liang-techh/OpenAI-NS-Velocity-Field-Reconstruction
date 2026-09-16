from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_osc_covariance_increment_replay import (
    ExactEnsemble,
    bilinear_covariance,
    covariance_increment_audit,
    covariance_remainder_exponents,
    matrix_add,
    matrix_sub,
)


WEIGHTS = (F(1, 3), F(2, 3))
W0 = ExactEnsemble(WEIGHTS, ((F(1), F(2)), (F(-1), F(1, 2))))
E = ExactEnsemble(WEIGHTS, ((F(1, 5), F(-2, 7)), (F(3, 4), F(1, 6))))
L = ExactEnsemble(WEIGHTS, ((F(2, 3), F(-1, 4)), (F(1, 5), F(4, 3))))
R = ExactEnsemble(WEIGHTS, ((F(-1, 6), F(2, 5)), (F(2, 7), F(-3, 8))))
SIGMA = ((F(8, 45), F(-83, 180)), (F(-83, 180), F(5, 9)))


def test_exact_covariance_increment_and_quadratic_split():
    result = covariance_increment_audit(W0, E, L, R, SIGMA)
    assert result["delta_covariance"] == (
        (F(7037, 14700), F(18577, 25200)),
        (F(18577, 25200), F(248429, 151200)),
    )
    assert result["quadratic_remainder"] == matrix_add(
        result["linear_quadratic"],
        result["linear_curl_cross"],
        result["curl_quadratic"],
    )


def test_each_retained_remainder_is_detectable():
    result = covariance_increment_audit(W0, E, L, R, SIGMA)
    delta = result["delta_covariance"]
    without_error_linear = matrix_add(
        result["sigma"], result["actual_curl_remainder"], result["quadratic_remainder"]
    )
    without_actual_curl = matrix_add(
        result["sigma"], result["error_linear_remainder"], result["quadratic_remainder"]
    )
    without_quadratic = matrix_add(
        result["sigma"], result["error_linear_remainder"], result["actual_curl_remainder"]
    )
    assert matrix_sub(delta, without_error_linear) == (
        (F(13, 45), F(767, 1260)),
        (F(767, 1260), F(65, 189)),
    )
    assert matrix_sub(delta, without_actual_curl) == (
        (F(-8, 35), F(6407, 25200)),
        (F(6407, 25200), F(13, 105)),
    )
    assert matrix_sub(delta, without_quadratic) == (
        (F(1179, 4900), F(169, 504)),
        (F(169, 504), F(13387, 21600)),
    )


def test_actual_u_must_be_retained_in_bilinear_curl_remainder():
    result = covariance_increment_audit(W0, E, L, R, SIGMA)
    frozen_only = bilinear_covariance(W0, R)
    defect = matrix_sub(result["actual_curl_remainder"], frozen_only)
    assert defect == bilinear_covariance(E, R)
    assert defect != ((F(0), F(0)), (F(0), F(0)))


def test_signed_inverse_premise_fails_closed_under_exact_mutation():
    mutated = ExactEnsemble(
        WEIGHTS,
        ((F(2, 3) + F(1, 2**40), F(-1, 4)), (F(1, 5), F(4, 3))),
    )
    defect = matrix_sub(bilinear_covariance(W0, mutated), SIGMA)
    assert defect == (
        (F(1, 1649267441664), F(1, 1649267441664)),
        (F(1, 1649267441664), F(0)),
    )
    with pytest.raises(ValueError, match="signed-inverse premise"):
        covariance_increment_audit(W0, E, mutated, R, SIGMA)


def test_displayed_remainder_exponents_are_exact():
    assert covariance_remainder_exponents(F(9, 10), F(3, 4), F(1, 100000)) == {
        "B(E,L)": F(23, 20),
        "B(U,R)": F(139999, 100000),
        "C(L)": F(4, 5),
        "B(L,R)": F(129999, 100000),
        "C(R)": F(89999, 50000),
    }


def test_approximate_theorem_inputs_are_rejected():
    with pytest.raises(TypeError):
        covariance_remainder_exponents(F(9, 10), 0.75, F(1, 100000))
    with pytest.raises(TypeError):
        ExactEnsemble((F(1, 2), 0.5), ((F(1), F(0)), (F(0), F(1))))


def test_weight_and_parameter_boundaries_fail_closed():
    with pytest.raises(ValueError, match="sum exactly to one"):
        ExactEnsemble((F(1, 3), F(1, 3)), ((F(1),), (F(2),)))
    with pytest.raises(ValueError, match="positive"):
        covariance_remainder_exponents(F(9, 10), F(3, 4), F(0))
