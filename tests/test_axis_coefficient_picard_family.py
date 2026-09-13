from decimal import Decimal, localcontext

import pytest

from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_mixed_scale import (
    MixedScaleCoefficient,
)
from openai_ns_reconstruction.axis_coefficient_picard_family import (
    formal_axis_picard_family_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_remainder import (
    wide_first_picard_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def states():
    x1 = actual_schedule_wide_first_picard_state(_data(), 0.05)
    solver = formal_axis_coefficient_solver(x1)
    family = formal_axis_picard_family_state(solver)
    return x1, solver, family, wide_first_picard_remainder_state(x1)


def _map_from_reference(solver, n: int, m: int, eta: float, angular: bool) -> MixedScaleCoefficient:
    phi, u = solver.x1.reference.jet_pair(n, m, eta)
    value = phi if angular else u
    return MixedScaleCoefficient.channel(0, 0, Decimal.from_float(value))


def _map_from_x1(x1, n: int, m: int, eta: float, angular: bool) -> MixedScaleCoefficient:
    phi, u = x1.jet_pair(n, m, eta)
    typed = phi if angular else u
    channels = {
        (0, 0): typed.reference,
        (0, 1): typed.inverse_lambda_numerator,
        (0, 2): typed.inverse_lambda_squared_numerator,
    }
    if not angular:
        pressure = x1.remainder.axial.pressure_normalized_factor(n, m, eta)
        with localcontext() as context:
            context.prec = 96
            channels[(2, 1)] = +(pressure / Decimal(2))
    return MixedScaleCoefficient(channels)


def _map_from_x2(jet) -> MixedScaleCoefficient:
    channels = {
        (0, power): value
        for power, value in enumerate(jet.ordinary)
        if value != 0
    }
    channels.update(
        {
            (2, power): value
            for power, value in enumerate(jet.pressure_linear)
            if value != 0
        }
    )
    if jet.pressure_square_inverse_lambda_fourth != 0:
        channels[(4, 4)] = jet.pressure_square_inverse_lambda_fourth
    return MixedScaleCoefficient(channels)


def _assert_sparse_close(
    actual: MixedScaleCoefficient,
    expected: MixedScaleCoefficient,
) -> None:
    with localcontext() as context:
        context.prec = 96
        for key in set(actual.support) | set(expected.support):
            left = actual.channels.get(key, Decimal(0))
            right = expected.channels.get(key, Decimal(0))
            scale = max(abs(left), abs(right), Decimal(1))
            assert abs(left - right) <= Decimal("5e-12") * scale


def test_picard_family_anchor_and_metadata(states) -> None:
    x1, solver, family, _ = states
    eta = 0.17
    rows = family.jet_prefix(0, 2, 0, eta)
    for n, (angular, axial) in enumerate(rows):
        _assert_sparse_close(angular, _map_from_reference(solver, n, 0, eta, True))
        _assert_sparse_close(axial, _map_from_reference(solver, n, 0, eta, False))
    assert family.Lambda == solver.Lambda == x1.Lambda
    assert family.finite_picard_materialized is True
    assert family.formal_coefficients_materialized is True
    assert family.fixed_point_materialized is False
    assert family.fixed_point_convergence_certified is False
    assert family.global_axis_norm_certified is False
    assert family.truncation_certified is False
    assert family.paper_exact is False


def test_picard_iterations_match_landed_x1_and_x2_adapters(states) -> None:
    x1, _, family, remainder = states
    eta = 0.17
    for m in (0, 1):
        first_rows = family.jet_prefix(1, 2, m, eta)
        second_rows = family.jet_prefix(2, 2, m, eta)
        for n in range(3):
            expected_first = (
                _map_from_x1(x1, n, m, eta, True),
                _map_from_x1(x1, n, m, eta, False),
            )
            expected_second_jets = remainder.second_picard_jet_pair(n, m, eta)
            expected_second = tuple(_map_from_x2(jet) for jet in expected_second_jets)
            _assert_sparse_close(first_rows[n][0], expected_first[0])
            _assert_sparse_close(first_rows[n][1], expected_first[1])
            _assert_sparse_close(second_rows[n][0], expected_second[0])
            _assert_sparse_close(second_rows[n][1], expected_second[1])

        # Include the next row once so the adapter comparison covers n = 3
        # without turning every low-order check into a high-row scan.
        if m == 0:
            expected_second_jets = remainder.second_picard_jet_pair(3, m, eta)
            expected_second = tuple(_map_from_x2(jet) for jet in expected_second_jets)
            actual_second = family.jet_pair(2, 3, m, eta)
            _assert_sparse_close(actual_second[0], expected_second[0])
            _assert_sparse_close(actual_second[1], expected_second[1])


def test_picard_rows_stabilize_against_formal_solver(states) -> None:
    _, solver, family, _ = states
    eta = 0.17
    for m in (0, 1):
        first_rows = family.jet_prefix(1, 2, m, eta)
        second_rows = family.jet_prefix(2, 2, m, eta)
        third_rows = family.jet_prefix(3, 2, m, eta)
        for n in (0, 1, 2):
            assert second_rows[n] == third_rows[n]
        for n in (0, 1):
            expected = solver.jet_pair(n, m, eta)
            _assert_sparse_close(first_rows[n][0], expected[0])
            _assert_sparse_close(first_rows[n][1], expected[1])
        for n in (0, 1, 2):
            expected = solver.jet_pair(n, m, eta)
            _assert_sparse_close(second_rows[n][0], expected[0])
            _assert_sparse_close(second_rows[n][1], expected[1])
