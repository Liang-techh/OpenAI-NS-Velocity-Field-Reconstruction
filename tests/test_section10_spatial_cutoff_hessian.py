import numpy as np

from openai_ns_reconstruction.spatial_localization import (
    section10_spatial_cutoff_gradient,
    section10_spatial_cutoff_hessian,
)


def _gradient_difference_hessian(point: np.ndarray, h: float = 1.0e-6) -> np.ndarray:
    """Independent central difference of the already-existing analytic gradient."""
    columns = []
    for axis in range(3):
        step = np.zeros(3, dtype=float)
        step[axis] = h
        plus = section10_spatial_cutoff_gradient(*(point + step))
        minus = section10_spatial_cutoff_gradient(*(point - step))
        columns.append((plus - minus) / (2.0 * h))
    return np.column_stack(columns)


def test_section10_spatial_cutoff_hessian_matches_gradient_difference_in_collar():
    # Both arguments 16 r^2 and 4 z lie strictly inside the transition collar
    # (1/2, 1), so this exercises radial, axial, and mixed second derivatives.
    point = np.array([0.20, 0.02, 0.16], dtype=float)
    analytic = section10_spatial_cutoff_hessian(*point)
    independent = _gradient_difference_hessian(point)

    np.testing.assert_allclose(analytic, independent, rtol=2.0e-5, atol=2.0e-5)
    np.testing.assert_allclose(analytic, analytic.T, rtol=0.0, atol=0.0)


def test_section10_spatial_cutoff_hessian_is_zero_on_plateau_and_outside_support():
    np.testing.assert_array_equal(
        section10_spatial_cutoff_hessian(0.0, 0.0, 0.0), np.zeros((3, 3))
    )
    np.testing.assert_array_equal(
        section10_spatial_cutoff_hessian(0.30, 0.0, 0.0), np.zeros((3, 3))
    )
    np.testing.assert_array_equal(
        section10_spatial_cutoff_hessian(0.0, 0.0, 0.30), np.zeros((3, 3))
    )


def test_section10_spatial_cutoff_hessian_rejects_nonfinite_coordinates():
    for bad in (float("nan"), float("inf"), -float("inf")):
        try:
            section10_spatial_cutoff_hessian(bad, 0.0, 0.0)
        except ValueError:
            pass
        else:
            raise AssertionError("nonfinite coordinate should fail closed")
