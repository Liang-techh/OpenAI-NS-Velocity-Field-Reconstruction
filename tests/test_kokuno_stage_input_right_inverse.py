from decimal import Decimal
from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_stage_input_right_inverse import (
    assemble_squared_partition,
    covariance_response_from_deltas,
    invert_matrix2,
    signed_differential_inverse,
)


H = ((F(2), F(1)), (F(1), F(3)))
SIGMA = (F(7, 10), F(-2, 5))
EPSILON = F(1, 16)
AMPLITUDES = (F(3, 2), F(5, 4))


def test_signed_differential_inverse_recovers_covariance_exactly():
    receipt = signed_differential_inverse(H, SIGMA, EPSILON, AMPLITUDES)
    assert receipt.d_sigma == (F(8), F(-24, 5))
    assert receipt.delta_amplitudes == (F(8, 3), F(-48, 25))
    assert receipt.component_variation == (F(1, 2), F(-3, 10))
    assert receipt.covariance_response == SIGMA
    assert receipt.residual == (F(0), F(0))
    assert receipt.verified


def test_exact_inverse_is_not_assumed_and_singular_map_fails_closed():
    assert invert_matrix2(H) == ((F(3, 5), F(-1, 5)), (F(-1, 5), F(2, 5)))
    with pytest.raises(ValueError, match="invertible"):
        signed_differential_inverse(((1, 2), (2, 4)), SIGMA, EPSILON, AMPLITUDES)


def test_amplitude_delta_mutation_remains_exactly_visible():
    receipt = signed_differential_inverse(H, SIGMA, EPSILON, AMPLITUDES)
    mutated = (receipt.delta_amplitudes[0] + F(1, 2**40), receipt.delta_amplitudes[1])
    response = covariance_response_from_deltas(H, EPSILON, AMPLITUDES, mutated)
    residual = (response[0] - SIGMA[0], response[1] - SIGMA[1])
    assert residual == (F(3, 2**43), F(3, 2**44))
    assert residual != (0, 0)


def test_squared_partition_assembles_global_covariance_exactly():
    assembled = assemble_squared_partition(SIGMA, (F(3, 5), F(4, 5)))
    assert assembled.square_sum == 1
    assert assembled.assembled == SIGMA
    assert assembled.residual == (0, 0)
    assert assembled.verified


def test_partition_weight_mutation_fails_closed_at_zero_tolerance():
    assembled = assemble_squared_partition(SIGMA, (F(3, 5) + F(1, 2**40), F(4, 5)))
    assert assembled.square_sum != 1
    assert assembled.residual != (0, 0)
    assert not assembled.verified


def test_nonpositive_inputs_fail_closed():
    with pytest.raises(ValueError, match="epsilon"):
        signed_differential_inverse(H, SIGMA, 0, AMPLITUDES)
    with pytest.raises(ValueError, match="amplitudes"):
        signed_differential_inverse(H, SIGMA, EPSILON, (0, 1))
    with pytest.raises(ValueError, match="partition"):
        assemble_squared_partition(SIGMA, ())


def test_approximate_theorem_inputs_are_rejected():
    with pytest.raises(TypeError):
        signed_differential_inverse(H, SIGMA, 0.0625, AMPLITUDES)
    with pytest.raises(TypeError):
        signed_differential_inverse(H, SIGMA, Decimal("0.0625"), AMPLITUDES)
    with pytest.raises(TypeError):
        assemble_squared_partition(SIGMA, (0.6, F(4, 5)))
