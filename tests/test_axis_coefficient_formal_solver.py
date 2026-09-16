from decimal import Decimal, localcontext

import pytest

from openai_ns_reconstruction.axis_coefficient_mixed_scale import (
    MixedScaleCoefficient,
    mixed_scale_product,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_remainder import (
    wide_first_picard_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)


def _polynomial_family(q: int, p: int, constant: str, slope: str):
    constant = Decimal(constant)
    slope = Decimal(slope)

    def family(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        if n != 0:
            return MixedScaleCoefficient.zero()
        eta_decimal = Decimal(str(eta))
        if m == 0:
            value = constant + slope * eta_decimal
        elif m == 1:
            value = slope
        else:
            value = Decimal(0)
        return MixedScaleCoefficient.channel(q, p, value)

    return family


def test_product_uses_eta_binomial_and_keeps_a4_channel() -> None:
    left = _polynomial_family(2, 1, "2", "3")
    right = _polynomial_family(2, 2, "5", "7")
    product = mixed_scale_product(left, right)(0, 1, 0.25)
    assert product.support == ((4, 3),)
    # (2 + 3 eta) * 7 + 3 * (5 + 7 eta) at eta = 1/4.
    assert product[(4, 3)] == Decimal("39.50")


def _second_picard_map(jet) -> MixedScaleCoefficient:
    channels = {}
    channels.update(
        ((0, power), value)
        for power, value in enumerate(jet.ordinary)
        if value != 0
    )
    channels.update(
        ((2, power), value)
        for power, value in enumerate(jet.pressure_linear)
        if value != 0
    )
    if jet.pressure_square_inverse_lambda_fourth != 0:
        channels[(4, 4)] = jet.pressure_square_inverse_lambda_fourth
    return MixedScaleCoefficient(channels)


def _assert_sparse_close(actual: MixedScaleCoefficient, expected: MixedScaleCoefficient) -> None:
    keys = set(actual.support) | set(expected.support)
    with localcontext() as ctx:
        ctx.prec = 96
        for key in keys:
            left = actual.channels.get(key, Decimal(0))
            right = expected.channels.get(key, Decimal(0))
            scale = max(abs(left), abs(right), Decimal(1))
            assert abs(left - right) <= Decimal("5e-12") * scale


def test_formal_solver_row_zero_keeps_reference_and_truth_boundary(x1) -> None:
    state = formal_axis_coefficient_solver(x1)
    angular, axial = state.jet_pair(0, 0, 0.17)
    assert angular.channels.get((0, 0), Decimal(0)) == Decimal.from_float(
        x1.reference.phi.jet(0, 0, 0.17)
    )
    assert axial.channels.get((0, 0), Decimal(0)) == Decimal.from_float(
        x1.reference.u.jet(0, 0, 0.17)
    )
    assert state.fixed_point_materialized is False
    assert state.global_axis_norm_certified is False
    assert state.paper_exact is False


def test_formal_solver_rows_match_second_picard_channelwise(x1) -> None:
    state = formal_axis_coefficient_solver(x1)
    remainder = wide_first_picard_remainder_state(x1)
    for n in (1, 2):
        actual_angular, actual_axial = state.jet_pair(n, 0, 0.17)
        expected_angular, expected_axial = remainder.second_picard_jet_pair(n, 0, 0.17)
        _assert_sparse_close(actual_angular, _second_picard_map(expected_angular))
        _assert_sparse_close(actual_axial, _second_picard_map(expected_axial))
