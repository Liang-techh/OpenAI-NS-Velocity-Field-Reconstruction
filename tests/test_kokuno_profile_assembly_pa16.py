from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_profile_assembly_pa16 import (
    audit_pa16_power_blocks,
)


def test_pa16_exact_fixture_is_invertible() -> None:
    audit = audit_pa16_power_blocks(
        u_points_y=(Fraction(1, 2), Fraction(3, 4)),
        e_points_y=(Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)),
    )
    assert audit.u_determinant == Fraction(665, 4096)
    assert audit.e_determinant == Fraction(137274860048095, 5015306502144)
    assert audit.u_sign == 1
    assert audit.e_sign == 1
    assert audit.invertible


def test_pa16_exact_small_separation_remains_nonzero() -> None:
    delta = Fraction(1, 2**40)
    audit = audit_pa16_power_blocks(
        u_points_y=(Fraction(1, 2), Fraction(1, 2) + delta),
        e_points_y=(Fraction(1, 2), Fraction(1, 2) + delta, Fraction(3, 4)),
    )
    assert audit.u_determinant > 0
    assert audit.e_determinant > 0


def test_pa16_collision_fails_closed() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_pa16_power_blocks(
            u_points_y=(Fraction(1, 2), Fraction(1, 2)),
            e_points_y=(Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)),
        )


def test_pa16_e_block_collision_fails_closed() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_pa16_power_blocks(
            u_points_y=(Fraction(1, 2), Fraction(3, 4)),
            e_points_y=(Fraction(1, 2), Fraction(2, 3), Fraction(2, 3)),
        )


def test_pa16_order_reversal_fails_closed() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_pa16_power_blocks(
            u_points_y=(Fraction(3, 4), Fraction(1, 2)),
            e_points_y=(Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)),
        )


def test_pa16_nonpositive_support_fails_closed() -> None:
    with pytest.raises(ValueError, match="positive"):
        audit_pa16_power_blocks(
            u_points_y=(Fraction(0), Fraction(1, 2)),
            e_points_y=(Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)),
        )


def test_pa16_float_theorem_data_fails_closed() -> None:
    with pytest.raises(TypeError, match="Fraction"):
        audit_pa16_power_blocks(
            u_points_y=(Fraction(1, 2), 0.75),
            e_points_y=(Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)),
        )
