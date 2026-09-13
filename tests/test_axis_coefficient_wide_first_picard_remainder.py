from decimal import Decimal, localcontext

import pytest

from openai_ns_reconstruction.axis_coefficient_data import (
    actual_schedule_axis_coefficient_data,
)
from openai_ns_reconstruction.axis_coefficient_natural_resolvent import (
    actual_schedule_natural_resolvent,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_remainder import (
    wide_first_picard_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def remainder(x1):
    return wide_first_picard_remainder_state(x1)


def _assert_remainder_shape(jet) -> None:
    assert len(jet.ordinary) == 6
    assert len(jet.pressure_linear) == 5
    assert all(isinstance(value, Decimal) for value in jet.ordinary)
    assert all(isinstance(value, Decimal) for value in jet.pressure_linear)
    assert isinstance(jet.pressure_square_inverse_lambda_cubed, Decimal)


def _assert_second_shape(jet) -> None:
    assert len(jet.ordinary) == 7
    assert len(jet.pressure_linear) == 6
    assert all(isinstance(value, Decimal) for value in jet.ordinary)
    assert all(isinstance(value, Decimal) for value in jet.pressure_linear)
    assert isinstance(jet.pressure_square_inverse_lambda_fourth, Decimal)


def _assert_decimal_close(actual: Decimal, expected: Decimal) -> None:
    """Allow the x0 binary64 oracle's conversion error at a narrow bound."""

    with localcontext() as ctx:
        ctx.prec = PRECISION
        scale = max(abs(actual), abs(expected), Decimal(1))
        assert abs(actual - expected) <= Decimal("5e-12") * scale


def test_remainder_metadata_and_exact_zero_row(remainder, x1) -> None:
    assert remainder.x1 is x1
    assert remainder.epsilon == x1.epsilon
    assert remainder.Lambda == x1.Lambda
    assert remainder.natural_remainder_x1_materialized is True
    assert remainder.fixed_point_materialized is False
    assert remainder.fixed_point_convergence_certified is False
    assert remainder.paper_exact is False
    assert remainder.global_axis_norm_certified is False

    angular, axial = remainder.jet_pair(0, 0, 0.0)
    for jet in (angular, axial):
        _assert_remainder_shape(jet)
        assert jet.ordinary == (Decimal(0),) * 6
        assert jet.pressure_linear == (Decimal(0),) * 5
        assert jet.pressure_square_inverse_lambda_cubed == Decimal(0)
        assert jet.Lambda == x1.Lambda


def test_power_zero_channels_match_landed_x0_remainder(remainder, x1) -> None:
    pressure_seen = False
    for n in (1, 2):
        angular, axial = remainder.jet_pair(n, 0, 0.0)
        expected_angular = Decimal.from_float(
            x1.remainder.angular.ordinary_base.jet(n, 0, 0.0)
        )
        expected_axial = Decimal.from_float(
            x1.remainder.axial.ordinary_base.jet(n, 0, 0.0)
        )
        _assert_decimal_close(angular.ordinary[0], expected_angular)
        _assert_decimal_close(axial.ordinary[0], expected_axial)

        expected_pressure = x1.remainder.axial.pressure_normalized_factor(n, 0, 0.0)
        _assert_decimal_close(axial.pressure_linear[0], expected_pressure)
        pressure_seen |= expected_pressure != 0

    # The first nonzero actual pressure channel is already visible at this
    # low row; this guards against an aggregate that silently drops pressure.
    assert pressure_seen


def test_second_picard_keeps_channel_arrays_under_outer_scale(remainder, x1) -> None:
    n, m, eta = 2, 0, 0.0
    source_angular, source_axial = remainder.jet_pair(n, m, eta)
    angular, axial = remainder.second_picard_jet_pair(n, m, eta)
    _assert_second_shape(angular)
    _assert_second_shape(axial)

    with localcontext() as ctx:
        ctx.prec = PRECISION
        half = Decimal(1) / Decimal(2)
        expected_angular_ordinary = (
            Decimal.from_float(x1.reference.phi.jet(n, m, eta)),
        ) + tuple(+(value * half) for value in source_angular.ordinary)
        expected_axial_ordinary = (
            Decimal.from_float(x1.reference.u.jet(n, m, eta)),
        ) + tuple(+(value * half) for value in source_axial.ordinary)
        expected_angular_pressure = (Decimal(0),) + tuple(
            +(value * half) for value in source_angular.pressure_linear
        )
        expected_axial_pressure = (Decimal(0),) + tuple(
            +(value * half) for value in source_axial.pressure_linear
        )
        assert angular.ordinary == expected_angular_ordinary
        assert axial.ordinary == expected_axial_ordinary
        assert angular.pressure_linear == expected_angular_pressure
        assert axial.pressure_linear == expected_axial_pressure
        assert angular.pressure_square_inverse_lambda_fourth == +(
            source_angular.pressure_square_inverse_lambda_cubed * half
        )
        assert axial.pressure_square_inverse_lambda_fourth == +(
            source_axial.pressure_square_inverse_lambda_cubed * half
        )

    assert angular.Lambda == axial.Lambda == x1.Lambda
    expected_log = x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
    assert angular.amplitude_log == axial.amplitude_log == expected_log


def test_natural_resolvent_triangle_once(x1) -> None:
    """Check linearity on one small actual coefficient triangle."""

    reference = x1.reference
    data = actual_schedule_axis_coefficient_data(reference)
    resolvent = actual_schedule_natural_resolvent(reference)
    left, right = data.one, data.eta
    summed = AxisCoefficientJetState(
        epsilon=reference.epsilon,
        origin="test-only actual resolvent triangle sum",
        _jet_provider=lambda n, m, eta: left.jet(n, m, eta) + right.jet(n, m, eta),
    )
    resolved_sum = resolvent(summed)
    separate_left = resolvent(left)
    separate_right = resolvent(right)
    n, m, eta = 2, 1, 0.17
    assert resolved_sum.jet(n, m, eta) == pytest.approx(
        separate_left.jet(n, m, eta) + separate_right.jet(n, m, eta),
        rel=2e-9,
        abs=5e-13,
    )
