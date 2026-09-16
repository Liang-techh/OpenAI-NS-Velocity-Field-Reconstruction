from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_mean_mc30_gain_floor import (
    SOURCE_KAPPA_S,
    replay_exposed_mc30_gain_floor,
    source_exposed_mc30_exponents,
    verify_exposed_mc30_gain_floor,
)


def test_boundary_alpha_exact_source_margins():
    out = replay_exposed_mc30_gain_floor(Fraction(9, 10))
    assert out.target_exponent == Fraction(89999, 50000)
    assert dict(out.term_exponents) == {
        "rg_time_epsilon_dT_B": Fraction(29, 10),
        "rg_radial_2bB": Fraction(289999, 100000),
        "rg_radial_2betaB": Fraction(379999, 100000),
        "rg_radial_B2": Fraction(379999, 100000),
        "rg_axial_bd": Fraction(29, 10),
        "mc29_product_alpha_plus_09": Fraction(9, 5),
        "mc29_product_2alpha": Fraction(9, 5),
    }
    assert dict(out.margins)["mc29_product_alpha_plus_09"] == Fraction(1, 50000)
    assert dict(out.margins)["mc29_product_2alpha"] == Fraction(1, 50000)
    assert out.minimum_margin == Fraction(1, 50000)


def test_larger_alpha_still_meets_floor_exactly():
    out = replay_exposed_mc30_gain_floor(Fraction(7, 5))
    assert all(margin >= 0 for _, margin in out.margins)
    assert out.minimum_margin == 2 * SOURCE_KAPPA_S


def test_one_extra_kappa_loss_in_product_row_fails_closed():
    alpha = Fraction(9, 10)
    rows = source_exposed_mc30_exponents(alpha)
    target = alpha + Fraction(9, 10) - 2 * SOURCE_KAPPA_S
    rows["mc29_product_alpha_plus_09"] = target - SOURCE_KAPPA_S
    with pytest.raises(ValueError, match="target exponent"):
        verify_exposed_mc30_gain_floor(alpha, rows)


def test_missing_or_extra_rows_fail_closed():
    alpha = Fraction(9, 10)
    rows = source_exposed_mc30_exponents(alpha)
    rows.pop("rg_axial_bd")
    with pytest.raises(ValueError, match="keys mismatch"):
        verify_exposed_mc30_gain_floor(alpha, rows)

    rows = source_exposed_mc30_exponents(alpha)
    rows["invented"] = Fraction(99)
    with pytest.raises(ValueError, match="keys mismatch"):
        verify_exposed_mc30_gain_floor(alpha, rows)


def test_source_hypothesis_alpha_boundary_is_fail_closed():
    with pytest.raises(ValueError, match="alpha >= 9/10"):
        replay_exposed_mc30_gain_floor(Fraction(9, 10) - Fraction(1, 10**12))


def test_source_kappa_is_pinned_exactly():
    with pytest.raises(ValueError, match="published source value"):
        replay_exposed_mc30_gain_floor(Fraction(9, 10), Fraction(1, 99999))


def test_approximate_theorem_inputs_rejected():
    with pytest.raises(TypeError):
        replay_exposed_mc30_gain_floor(0.9)  # type: ignore[arg-type]
    rows = source_exposed_mc30_exponents(Fraction(9, 10))
    rows["rg_axial_bd"] = 2.9  # type: ignore[assignment]
    with pytest.raises(TypeError):
        verify_exposed_mc30_gain_floor(Fraction(9, 10), rows)
