import math

import pytest

from openai_ns_reconstruction.background_picard_bounds import (
    eq_5_8_log_majorant,
    eq_5_8_majorant,
    find_picard_truncation,
    picard_parameter_derivative_count,
    picard_tail_certificate,
)


def test_parameter_derivative_count_is_exact_ceil_k_over_two():
    assert [picard_parameter_derivative_count(k) for k in range(8)] == [
        0,
        1,
        1,
        2,
        2,
        3,
        3,
        4,
    ]


def test_eq_5_8_majorant_matches_displayed_formula():
    # k=4 gives p_k=2.  Here C_n*a=1 and Delta=1/4, so the
    # Cauchy factor is 8^2 and the simplex factor is 1/5!.
    got = eq_5_8_majorant(4, C_n=2.0, radial_extent=0.5, rho=1.0, rho_prime=0.75)
    expected = 64.0 / math.factorial(5)
    assert got >= expected
    assert got == pytest.approx(expected, rel=2e-15)
    assert math.log(got) == pytest.approx(
        eq_5_8_log_majorant(4, 2.0, 0.5, 1.0, 0.75), rel=2e-15
    )


def test_eq_5_8_majorant_dominates_independent_scalar_picard_specialization():
    """Cross-check the simplex factor against an exactly solvable K iteration.

    Take one c_i=0 coordinate, A1=0, A0=lambda and f=1.  Then repeated radial
    integration gives the exact k-th term

        K^k Gf(xi) = lambda^k xi^(k+1)/(k+1)!.

    This oracle is derived without calling the production recurrence or bound.
    """

    lam = 0.75
    xi = 0.6
    radial_extent = 0.8
    # C_n=1 bounds both |lambda| and |f| in this scalar specialization.
    for k in range(12):
        exact_term = lam**k * xi ** (k + 1) / math.factorial(k + 1)
        bound = eq_5_8_majorant(k, 1.0, radial_extent, 1.0, 0.5)
        assert exact_term <= bound


def test_two_step_tail_certificate_dominates_direct_majorant_partial_tail():
    cert = picard_tail_certificate(
        8, C_n=2.0, radial_extent=0.5, rho=1.0, rho_prime=0.5
    )
    # Independently sum many terms of the displayed Eq. (5.8) expression.  The
    # production tail certificate instead uses the analytic two-step ratio.
    direct = 0.0
    for k in range(8, 160):
        p = (k + 1) // 2
        direct += (
            1.0 ** (k + 1)
            / math.factorial(k + 1)
            * max(1.0, p / 0.5) ** p
        )
    assert direct <= cert.tail_upper_bound
    assert 0.0 < cert.two_step_ratio_upper < 1.0


def test_find_picard_truncation_is_bound_driven_and_first_admissible_for_fixture():
    cert = find_picard_truncation(
        1.0e-8,
        C_n=2.0,
        radial_extent=0.5,
        rho=1.0,
        rho_prime=0.5,
        max_start_order=100,
    )
    assert cert.start_order == 24
    assert cert.certifies_tolerance(1.0e-8)

    previous = picard_tail_certificate(23, 2.0, 0.5, 1.0, 0.5)
    assert not previous.certifies_tolerance(1.0e-8)


def test_picard_tail_certificate_fails_closed_before_ratio_regime():
    # For these values the conservative two-step ratio at k=1 is exactly 1.
    with pytest.raises(ValueError, match="two-step ratio"):
        picard_tail_certificate(1, 2.0, 0.5, 1.0, 0.5)

    # A wide strip loss can also put p_k before the power branch used in the
    # tail derivation.
    with pytest.raises(ValueError, match="Cauchy-factor power regime"):
        picard_tail_certificate(2, 1.0, 0.5, 5.0, 1.0)


@pytest.mark.parametrize(
    "rho,rho_prime",
    [(1.0, 1.0), (1.0, 0.0), (1.0, -0.1), (float("nan"), 0.5)],
)
def test_eq_5_8_rejects_invalid_complex_strip(rho, rho_prime):
    with pytest.raises(ValueError):
        eq_5_8_majorant(0, 1.0, 1.0, rho, rho_prime)


def test_find_picard_truncation_refuses_unmet_search_cap():
    with pytest.raises(ValueError, match="no Eq. \(5.8\) Picard truncation"):
        find_picard_truncation(
            1.0e-30,
            C_n=2.0,
            radial_extent=0.5,
            rho=1.0,
            rho_prime=0.5,
            max_start_order=10,
        )
