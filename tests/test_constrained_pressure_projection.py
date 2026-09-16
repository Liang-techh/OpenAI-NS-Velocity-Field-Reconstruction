import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_projection import (
    bounded_pressure_projection,
)


def test_pressure_span_can_remove_only_declared_gradient_component():
    design = np.array([
        [1.0, 0.0],
        [0.0, 1.0],
        [0.0, 0.0],
    ])
    result = bounded_pressure_projection(
        np.array([-0.25, 0.5, 2.0]),
        design,
        lower_bounds=np.array([-1.0, -1.0]),
        upper_bounds=np.array([1.0, 1.0]),
        residual_scales=np.ones(3),
        parameter_labels=("px", "py"),
    )

    assert result.coefficients == pytest.approx((0.25, -0.5))
    assert result.residual_rms_before == pytest.approx(
        np.sqrt((0.25**2 + 0.5**2 + 2.0**2) / 3.0)
    )
    assert result.residual_rms_after == pytest.approx(2.0 / np.sqrt(3.0))
    assert result.design_rank == 2
    assert result.design_nullity == 0
    assert result.velocity_changed is False
    assert result.independent_validation_required is True


def test_fixed_forcing_is_immutable_and_not_fit():
    result = bounded_pressure_projection(
        np.array([3.0, -1.0]),
        np.array([[1.0], [0.0]]),
        lower_bounds=np.array([-1.0]),
        upper_bounds=np.array([1.0]),
        residual_scales=np.ones(2),
        fixed_forcing=np.array([1.0, -1.0]),
        parameter_labels=("p0",),
    )

    assert result.forcing_mode == "fixed"
    assert result.coefficients[0] == pytest.approx(-1.0, abs=1e-10)
    assert result.residual_rms_after == pytest.approx(1.0 / np.sqrt(2.0))
    assert result.active_lower_bounds == ("p0",)


def test_gauge_like_zero_column_is_reported_as_nullity():
    result = bounded_pressure_projection(
        np.array([-0.4, 0.3]),
        np.array([[1.0, 0.0], [0.0, 0.0]]),
        lower_bounds=-np.ones(2),
        upper_bounds=np.ones(2),
        residual_scales=np.ones(2),
        parameter_labels=("gradient_mode", "constant_gauge"),
    )

    assert result.design_rank == 1
    assert result.design_nullity == 1
    assert result.design_condition == pytest.approx(1.0)
    assert result.residual_rms_after == pytest.approx(abs(0.3) / np.sqrt(2.0))


def test_bounds_expose_pressure_capacity_limit():
    result = bounded_pressure_projection(
        np.array([-4.0]),
        np.array([[1.0]]),
        lower_bounds=np.array([-0.5]),
        upper_bounds=np.array([0.5]),
        residual_scales=np.ones(1),
        parameter_labels=("bounded_pressure",),
    )

    assert result.coefficients[0] == pytest.approx(0.5, abs=1e-10)
    assert result.residual_rms_after == pytest.approx(3.5, abs=1e-10)
    assert result.active_upper_bounds == ("bounded_pressure",)
    assert result.recoverable_fraction_rms == pytest.approx(0.125, abs=1e-10)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"pressure_gradient_design": np.ones((2, 2))}, "row count"),
        ({"residual_scales": np.array([1.0, 0.0, 1.0])}, "strictly positive"),
        (
            {"lower_bounds": np.array([0.0]), "upper_bounds": np.array([0.0])},
            "below",
        ),
        ({"fixed_forcing": np.array([1.0, np.nan, 0.0])}, "finite"),
    ],
)
def test_fail_closed_on_malformed_inputs(kwargs, match):
    common = dict(
        momentum_residual=np.ones(3),
        pressure_gradient_design=np.ones((3, 1)),
        lower_bounds=np.array([-1.0]),
        upper_bounds=np.array([1.0]),
        residual_scales=np.ones(3),
    )
    common.update(kwargs)
    with pytest.raises(ValueError, match=match):
        bounded_pressure_projection(**common)
