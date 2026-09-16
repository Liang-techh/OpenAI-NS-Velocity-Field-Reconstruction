from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_mean_vandermonde_replay import (
    SOURCE_ARCHIVE_PATH,
    SOURCE_ASSOCIATED_CHECKER,
    SOURCE_CHECKER_SHA256,
    SOURCE_COMPONENT,
    SOURCE_COMPONENT_SHA256,
    SOURCE_RESULT_SHA256,
    angular_interpolant,
    angular_source_determinant,
    angular_targets,
    axial_coefficients,
    axial_residuals,
    evaluate_quadratic,
    translated_moment,
)


def test_pinned_public_provenance() -> None:
    assert SOURCE_COMPONENT == "NS-mean_corrections"
    assert SOURCE_ARCHIVE_PATH == "proof_sources/mean_corrections/mean_corrections_body.tex"
    assert SOURCE_COMPONENT_SHA256 == "c494cebcaf476a01c5aae17d38f84dfc6b0092040c03b9aecbad1f5c3831ea29"
    assert SOURCE_ASSOCIATED_CHECKER == "proof_sources/mean_corrections/verify_exact.py"
    assert SOURCE_CHECKER_SHA256 == "5d54f54000e0cf03e70cfde1d298398bd29d091addf07887f70660af7c98a366"
    assert SOURCE_RESULT_SHA256 == "5b81365d43a7108b16d11eba8cdf1dc267266226a0d90ff221ff67827f37495f"


def test_mc23_translated_moment_scaling_exact() -> None:
    mu = F(7, 5)
    t = F(9, 4)
    assert translated_moment(mu, t, 0) == F(7, 5)
    assert translated_moment(mu, t, 1) == F(63, 20)
    assert translated_moment(mu, t, 2) == F(567, 80)


def test_mc25_angular_interpolation_and_determinant_exact() -> None:
    a = F(9, 8)
    mu0, mu1, mu2 = F(3, 2), F(5, 4), F(7, 6)
    t0, t1, t2 = F(2), F(3), F(5)
    b1, b2 = angular_targets(
        pressure_target=F(11, 7),
        axial_target=F(-13, 9),
        a=a,
        mu_p1=mu1,
        mu_p2=mu2,
    )
    coeffs = angular_interpolant(t0, t1, t2, b1, b2)
    assert evaluate_quadratic(coeffs, t0) == 0
    assert evaluate_quadratic(coeffs, t1) == b1
    assert evaluate_quadratic(coeffs, t2) == b2
    assert angular_source_determinant(
        a=a, mu_p0=mu0, mu_p1=mu1, mu_p2=mu2, t0=t0, t1=t1, t2=t2
    ) != 0


def test_mc25_coefficient_mutation_is_detected_with_zero_tolerance() -> None:
    coeffs = angular_interpolant(F(2), F(3), F(5), F(-4, 9), F(8, 13))
    mutated = (coeffs[0] + F(1, 2**40), coeffs[1], coeffs[2])
    assert evaluate_quadratic(mutated, F(2)) == F(1, 2**40)


def test_mc25_node_collision_fails_closed() -> None:
    with pytest.raises(ValueError, match="pairwise-distinct"):
        angular_interpolant(F(2), F(2), F(5), F(1), F(2))
    assert angular_source_determinant(
        a=F(1), mu_p0=F(1), mu_p1=F(1), mu_p2=F(1), t0=F(2), t1=F(2), t2=F(5)
    ) == 0


def test_mc26_axial_two_node_solve_exact() -> None:
    args = dict(t_b=F(4, 3), t_a=F(7, 5), j_theta=F(17, 11), a=F(9, 8), mu_weighted=F(6, 5))
    s0, s1 = axial_coefficients(**args)
    assert (s0, s1) == (F(6800, 297), F(-1700, 99))
    assert axial_residuals(
        s0=s0,
        s1=s1,
        t_b=args["t_b"],
        t_a=args["t_a"],
        j_theta=args["j_theta"],
        a=args["a"],
        mu_zero=F(5, 7),
        mu_weighted=args["mu_weighted"],
    ) == (0, 0)


def test_mc26_lambda_zero_style_collision_fails_closed() -> None:
    with pytest.raises(ValueError, match="lambda=0"):
        axial_coefficients(t_b=F(4, 3), t_a=F(4, 3), j_theta=F(1), a=F(1), mu_weighted=F(1))


def test_float_inputs_fail_closed() -> None:
    with pytest.raises(TypeError):
        translated_moment(F(1), 1.25, 1)  # type: ignore[arg-type]
