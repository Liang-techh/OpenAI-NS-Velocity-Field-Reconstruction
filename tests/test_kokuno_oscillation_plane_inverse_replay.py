from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_oscillation_plane_inverse_replay import (
    FULL_RECONSTRUCTION,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    MovingPlaneFixture,
    max_abs_entry,
    replay_plane_inverse,
)


def fixture() -> MovingPlaneFixture:
    return MovingPlaneFixture(
        s_a=F(2, 3),
        k_theta=F(3, 5),
        k_z=F(4, 5),
        tangential_norm=F(7, 5),
        gamma=F(5, 4),
    )


def test_exact_left_inverse_and_ambient_projector() -> None:
    result = replay_plane_inverse(fixture())
    assert result.passed
    assert max_abs_entry(result.left_inverse_residual) == 0
    assert max_abs_entry(result.ambient_formula_residual) == 0


def test_projector_has_source_range_and_kernel() -> None:
    result = replay_plane_inverse(fixture())
    assert max_abs_entry(result.projector_idempotence_residual) == 0
    assert max_abs_entry(result.normal_range_residual) == 0
    assert max_abs_entry(result.kernel_residual) == 0


def test_exact_mutation_of_reported_normal_fails_closed() -> None:
    data = fixture()
    eps = F(1, 2**40)
    bad_normal = (data.n_phi[0] + eps, data.n_phi[1], data.n_phi[2])
    result = replay_plane_inverse(data, projector_normal=bad_normal)
    assert max_abs_entry(result.ambient_formula_residual) == F(1, 7 * 2**38)
    assert not result.passed


def test_nonunit_tangential_direction_is_rejected() -> None:
    with pytest.raises(ValueError, match="unit tangential"):
        MovingPlaneFixture(F(2, 3), F(1, 2), F(1, 2), F(7, 5), F(5, 4))


def test_singular_coordinate_matrix_is_rejected() -> None:
    with pytest.raises(ValueError, match="invertible"):
        MovingPlaneFixture(F(2, 3), F(3, 5), F(4, 5), F(7, 5), F(0))


def test_non_fraction_inputs_are_rejected() -> None:
    with pytest.raises(TypeError, match="Fraction"):
        MovingPlaneFixture(2 / 3, F(3, 5), F(4, 5), F(7, 5), F(5, 4))  # type: ignore[arg-type]


def test_truth_boundary_stays_fail_closed() -> None:
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
