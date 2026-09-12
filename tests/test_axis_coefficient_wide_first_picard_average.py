from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_average import axis_coefficient_average
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_average import (
    ActualScheduleWideFirstPicardAverageState,
    MixedScaleFirstPicardAxialAverageCoefficientJet,
    wide_first_picard_average_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def states():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return x1, wide_first_picard_average_state(x1)


def test_average_state_is_bound_to_genuine_x1_and_keeps_truth_boundary(states) -> None:
    x1, average = states
    assert isinstance(average, ActualScheduleWideFirstPicardAverageState)
    assert average.x1 is x1
    assert average.Lambda == x1.Lambda
    assert average.epsilon == x1.epsilon
    assert average.mixed_scale_average_u1_materialized is True
    assert average.natural_remainder_x1_materialized is False
    assert average.fixed_point_materialized is False
    assert average.paper_exact is False


def test_average_scales_every_ordinary_x1_piece_by_exact_row_factor(states) -> None:
    x1, average = states
    n, m, eta = 3, 1, 0.04
    _, axial = x1.jet_pair(n, m, eta)
    averaged = average.jet(n, m, eta)

    assert isinstance(averaged, MixedScaleFirstPicardAxialAverageCoefficientJet)
    divisor = Decimal(n + 1)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        assert averaged.reference == +(axial.reference / divisor)
        assert averaged.inverse_lambda_numerator == +(
            axial.inverse_lambda_numerator / divisor
        )
        assert averaged.inverse_lambda_squared_numerator == +(
            axial.inverse_lambda_squared_numerator / divisor
        )
        first = +(
            axial.inverse_lambda_numerator
            / divisor
            / x1.Lambda
        )
        second = +(
            axial.inverse_lambda_squared_numerator
            / divisor
            / x1.Lambda
            / x1.Lambda
        )
    assert averaged.ordinary_correction_terms_decimal() == (first, second)


def test_average_reference_piece_cross_checks_landed_float_operator(states) -> None:
    x1, average = states
    reference_average = axis_coefficient_average(x1.reference.u)

    for n, m, eta in [(0, 0, 0.0), (1, 1, -0.03), (4, 0, 0.07)]:
        averaged = average.jet(n, m, eta)
        expected = reference_average.jet(n, m, eta)
        assert math.isclose(
            float(averaged.reference),
            expected,
            rel_tol=2e-15,
            abs_tol=2e-15,
        )


def test_average_pressure_keeps_signed_log_scale_and_exact_radial_divisor(states) -> None:
    x1, average = states
    eta = 0.0
    nonzero_row = next(
        (
            n
            for n in range(1, 12)
            if x1.jet_pair(n, 0, eta)[1].pressure_over_two_lambda.sign != 0
        ),
        None,
    )
    assert nonzero_row is not None

    _, axial = x1.jet_pair(nonzero_row, 0, eta)
    averaged = average.jet(nonzero_row, 0, eta)
    source = axial.pressure_over_two_lambda
    target = averaged.pressure_over_two_lambda_average

    assert source.sign != 0
    assert target.sign == source.sign
    assert target.log_scale == source.log_scale
    assert source.log_factor is not None
    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected_log_factor = +(
            source.log_factor - Decimal(nonzero_row + 1).ln()
        )
    assert target.log_factor == expected_log_factor

    # A nonzero theorem-scale pressure term may be too small for binary64, but
    # the lifted average must never silently replace it by zero.
    try:
        projected = target.to_binary64()
    except ArithmeticError as exc:
        assert "binary64" in str(exc)
    else:
        assert projected != 0.0
        assert math.isfinite(projected)


def test_row_zero_average_is_identity_on_exact_zero_pressure(states) -> None:
    x1, average = states
    _, axial = x1.jet_pair(0, 0, 0.0)
    averaged = average.jet(0, 0, 0.0)

    assert averaged.reference == axial.reference
    assert averaged.inverse_lambda_numerator == axial.inverse_lambda_numerator
    assert (
        averaged.inverse_lambda_squared_numerator
        == axial.inverse_lambda_squared_numerator
    )
    assert axial.pressure_over_two_lambda.sign == 0
    assert averaged.pressure_over_two_lambda_average.sign == 0


def test_average_fails_closed_for_wrong_state_or_negative_radial_degree(states) -> None:
    _, average = states
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        ActualScheduleWideFirstPicardAverageState(x1=object())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        average.jet(-1, 0, 0.0)


def test_average_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_average.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "average(u1)" in layer["capability"]
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert "fixed point" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
