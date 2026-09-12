import numpy as np
import pytest

from openai_ns_reconstruction.cutoffs import (
    standard_cutoff_second_derivative,
    standard_cutoff_third_derivative,
)
from openai_ns_reconstruction.section10_spatial_cutoff_viscous_adapter import (
    section10_spatial_cutoff_laplacian,
    section10_spatial_cutoff_laplacian_gradient,
    symmetric_smooth_bump_third_derivative,
)
from openai_ns_reconstruction.spatial_localization import (
    section10_spatial_cutoff_hessian,
)


def test_explicit_cutoff_third_derivative_matches_independent_second_difference():
    # Stay away from the flat joins so the finite-difference oracle samples one
    # smooth transition branch only.
    h = 1.0e-6
    for s in (0.61, 0.75, 0.89):
        independent = (
            standard_cutoff_second_derivative(s + h)
            - standard_cutoff_second_derivative(s - h)
        ) / (2.0 * h)
        assert standard_cutoff_third_derivative(s) == pytest.approx(
            independent, rel=2.0e-7, abs=2.0e-5
        )

    assert standard_cutoff_third_derivative(0.25) == 0.0
    assert standard_cutoff_third_derivative(1.25) == 0.0
    assert symmetric_smooth_bump_third_derivative(-0.75) == pytest.approx(
        -symmetric_smooth_bump_third_derivative(0.75)
    )


def _independent_gradient_of_hessian_trace(
    point: np.ndarray, h: float = 1.0e-6
) -> np.ndarray:
    result = np.empty(3, dtype=float)
    for axis in range(3):
        step = np.zeros(3, dtype=float)
        step[axis] = h
        plus = np.trace(section10_spatial_cutoff_hessian(*(point + step)))
        minus = np.trace(section10_spatial_cutoff_hessian(*(point - step)))
        result[axis] = (plus - minus) / (2.0 * h)
    return result


@pytest.mark.parametrize("z", [0.16, -0.16])
def test_grad_laplacian_matches_independent_hessian_trace_difference(z):
    # 16 r^2 and |4z| are both strictly inside (1/2, 1), exercising radial,
    # axial, and sign-sensitive third derivatives in the transition collar.
    point = np.array([0.20, 0.02, z], dtype=float)
    analytic = section10_spatial_cutoff_laplacian_gradient(*point)
    independent = _independent_gradient_of_hessian_trace(point)
    np.testing.assert_allclose(analytic, independent, rtol=2.0e-7, atol=2.0e-3)

    assert section10_spatial_cutoff_laplacian(*point) == pytest.approx(
        float(np.trace(section10_spatial_cutoff_hessian(*point))), rel=0.0, abs=0.0
    )


def test_cutoff_viscous_derivatives_are_exactly_zero_on_flat_regions():
    for point in ((0.0, 0.0, 0.0), (0.30, 0.0, 0.0), (0.0, 0.0, 0.30)):
        assert section10_spatial_cutoff_laplacian(*point) == 0.0
        np.testing.assert_array_equal(
            section10_spatial_cutoff_laplacian_gradient(*point), np.zeros(3)
        )


def test_cutoff_viscous_adapter_rejects_nonfinite_inputs():
    for bad in (float("nan"), float("inf"), -float("inf")):
        with pytest.raises(ValueError):
            section10_spatial_cutoff_laplacian_gradient(bad, 0.0, 0.0)
        with pytest.raises(ValueError):
            section10_spatial_cutoff_laplacian(0.0, bad, 0.0)
        with pytest.raises(ValueError):
            symmetric_smooth_bump_third_derivative(bad)
