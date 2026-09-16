from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_stage_input_covariance import (
    verify_stage_input_covariance,
)


def test_exact_stage_input_covariance_scaling_and_partition_square() -> None:
    receipt = verify_stage_input_covariance(
        h=Fraction(1, 8),
        A=Fraction(5, 8),
        partition_weights=(Fraction(3, 5), Fraction(4, 5)),
    )

    assert receipt.q_exponent == Fraction(-9, 8)
    assert receipt.Q_exponent == Fraction(0, 1)
    assert receipt.partition_square_sum == Fraction(1, 1)
    assert receipt.scaling_relation_ok is True
    assert receipt.partition_square_ok is True
    assert receipt.verified is True


def test_A_perturbation_is_not_tolerated() -> None:
    defect = Fraction(1, 2**40)
    receipt = verify_stage_input_covariance(
        h=Fraction(1, 8),
        A=Fraction(5, 8) + defect,
        partition_weights=(Fraction(3, 5), Fraction(4, 5)),
    )

    assert receipt.Q_exponent == -defect
    assert receipt.scaling_relation_ok is False
    assert receipt.verified is False


def test_partition_square_perturbation_is_not_tolerated() -> None:
    defect = Fraction(1, 2**40)
    receipt = verify_stage_input_covariance(
        h=Fraction(1, 8),
        A=Fraction(5, 8),
        partition_weights=(Fraction(3, 5) + defect, Fraction(4, 5)),
    )

    assert receipt.partition_square_sum != Fraction(1, 1)
    assert receipt.partition_square_ok is False
    assert receipt.verified is False


@pytest.mark.parametrize(
    ("h", "A", "weights"),
    [
        (0.125, Fraction(5, 8), (Fraction(3, 5), Fraction(4, 5))),
        (Fraction(1, 8), 0.625, (Fraction(3, 5), Fraction(4, 5))),
        (Fraction(1, 8), Fraction(5, 8), (0.6, Fraction(4, 5))),
    ],
)
def test_float_inputs_fail_closed(h, A, weights) -> None:
    with pytest.raises(TypeError):
        verify_stage_input_covariance(h=h, A=A, partition_weights=weights)


def test_empty_partition_does_not_verify() -> None:
    receipt = verify_stage_input_covariance(
        h=Fraction(1, 8),
        A=Fraction(5, 8),
        partition_weights=(),
    )

    assert receipt.partition_square_sum == Fraction(0, 1)
    assert receipt.partition_square_ok is False
    assert receipt.verified is False
