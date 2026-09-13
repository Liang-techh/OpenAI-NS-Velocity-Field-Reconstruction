from decimal import Decimal, localcontext
import math

import pytest

from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_mixed_scale import (
    MixedScaleCoefficient,
)
from openai_ns_reconstruction.axis_coefficient_profile_prefix import (
    evaluate_mixed_scale_radial_prefix,
    formal_axis_profile_prefix,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
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


def _row(scale: int) -> MixedScaleCoefficient:
    return MixedScaleCoefficient(
        {
            (0, 0): Decimal(2 * scale),
            (2, 1): Decimal(3 * scale),
            (4, 2): Decimal(5 * scale),
        }
    )


def _assert_sparse_close(
    actual: MixedScaleCoefficient,
    expected: MixedScaleCoefficient,
) -> None:
    with localcontext() as ctx:
        ctx.prec = 96
        for key in set(actual.support) | set(expected.support):
            left = actual.channels.get(key, Decimal(0))
            right = expected.channels.get(key, Decimal(0))
            scale = max(abs(left), abs(right), Decimal(1))
            assert abs(left - right) <= Decimal("1e-90") * scale


def test_hand_prefix_selects_exact_y0_rows_and_averages_before_derivative() -> None:
    rows = (_row(1), _row(2), _row(3))
    assert evaluate_mixed_scale_radial_prefix(rows, Decimal(0)) == rows[0]
    assert evaluate_mixed_scale_radial_prefix(rows, Decimal(0), radial_order=1) == rows[1]
    with localcontext() as ctx:
        ctx.prec = 96
        expected = rows[2].scale(+(Decimal(2) / Decimal(3)))
    assert evaluate_mixed_scale_radial_prefix(
        rows,
        Decimal(0),
        radial_order=2,
        average=True,
    ) == expected


def test_actual_prefix_matches_solver_rows_and_keeps_tail_unknown(x1) -> None:
    solver = formal_axis_coefficient_solver(x1)
    eta = 0.17
    for radial_order in (0, 1, 2):
        prefix = formal_axis_profile_prefix(
            solver,
            max_n=2,
            Y=Decimal(0),
            eta=eta,
            radial_order=radial_order,
            eta_order=0,
        )
        angular_row, axial_row = solver.jet_pair(radial_order, 0, eta)
        factor = Decimal(math.factorial(radial_order))
        assert prefix.angular == angular_row.scale(factor)
        assert prefix.axial == axial_row.scale(factor)
        with localcontext() as ctx:
            ctx.prec = 96
            expected_average = axial_row.scale(
                +(factor / Decimal(radial_order + 1))
            )
        _assert_sparse_close(prefix.axial_average, expected_average)

    prefix = formal_axis_profile_prefix(
        solver,
        max_n=2,
        Y=Decimal(1),
        eta=eta,
        radial_order=0,
        eta_order=0,
    )
    assert prefix.truncation_certified is False
    assert prefix.global_axis_norm_certified is False
    assert prefix.fixed_point_materialized is False
    assert prefix.paper_exact is False

    pressure_terms = [
        (key, term)
        for key, term in prefix.axial_terms_log().items()
        if key[0] > 0 and term.sign != 0
    ]
    assert pressure_terms
    _, pressure_term = pressure_terms[0]
    with localcontext() as ctx:
        ctx.prec = 96
        expected_log_scale = +(Decimal(2) * prefix.amplitude_log)
    assert pressure_term.log_scale == expected_log_scale
    try:
        projected = pressure_term.to_binary64()
    except ArithmeticError:
        pass
    else:
        assert projected != 0.0
