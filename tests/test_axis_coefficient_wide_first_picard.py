from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    ActualScheduleWideFirstPicardState,
    MixedScaleFirstPicardAngularCoefficientJet,
    MixedScaleFirstPicardAxialCoefficientJet,
    actual_schedule_wide_first_picard_state,
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


def test_first_picard_is_bound_to_complete_actual_remainder_and_scale(x1) -> None:
    assert isinstance(x1, ActualScheduleWideFirstPicardState)
    assert x1.reference is x1.remainder.reference
    assert x1.epsilon == x1.remainder.epsilon
    assert x1.Lambda == x1.remainder.Lambda
    assert x1.remainder.mixed_scale_natural_remainder_x0_complete is True
    assert x1.picard_x1_materialized is True
    assert x1.fixed_point_materialized is False
    assert x1.fixed_point_convergence_certified is False
    assert x1.paper_exact is False
    assert x1.global_axis_norm_certified is False


def test_angular_x1_uses_exact_outer_picard_scale_without_collapsing_into_x0(x1) -> None:
    n, m, eta = 2, 1, 0.07
    angular, _ = x1.jet_pair(n, m, eta)
    remainder = x1.remainder.angular.jet(n, m, eta)

    assert isinstance(angular, MixedScaleFirstPicardAngularCoefficientJet)
    assert angular.reference == Decimal.from_float(x1.reference.phi.jet(n, m, eta))
    assert angular.Lambda == remainder.Lambda == x1.Lambda
    with localcontext() as ctx:
        ctx.prec = PRECISION
        assert angular.inverse_lambda_numerator == +(remainder.ordinary_base / Decimal(2))
        assert angular.inverse_lambda_squared_numerator == +(
            remainder.inverse_lambda_numerator / Decimal(2)
        )
        expected_first = +(remainder.ordinary_base / Decimal(2) / x1.Lambda)
        expected_second = +(
            remainder.inverse_lambda_numerator
            / Decimal(2)
            / x1.Lambda
            / x1.Lambda
        )
    assert angular.correction_terms_decimal() == (expected_first, expected_second)


def test_axial_x1_preserves_ordinary_and_signed_log_pressure_corrections(x1) -> None:
    eta = 0.0
    nonzero_pressure_row = next(
        (
            n
            for n in range(12)
            if x1.remainder.axial.pressure_normalized_factor(n, 0, eta) != 0
        ),
        None,
    )
    assert nonzero_pressure_row is not None

    _, axial = x1.jet_pair(nonzero_pressure_row, 0, eta)
    remainder = x1.remainder.axial.jet(nonzero_pressure_row, 0, eta)
    assert isinstance(axial, MixedScaleFirstPicardAxialCoefficientJet)
    assert axial.reference == Decimal.from_float(
        x1.reference.u.jet(nonzero_pressure_row, 0, eta)
    )
    assert axial.Lambda == remainder.Lambda == x1.Lambda

    with localcontext() as ctx:
        ctx.prec = PRECISION
        assert axial.inverse_lambda_numerator == +(remainder.ordinary_base / Decimal(2))
        assert axial.inverse_lambda_squared_numerator == +(
            remainder.inverse_lambda_numerator / Decimal(2)
        )
        expected_first = +(remainder.ordinary_base / Decimal(2) / x1.Lambda)
        expected_second = +(
            remainder.inverse_lambda_numerator
            / Decimal(2)
            / x1.Lambda
            / x1.Lambda
        )
        log_divisor = +(Decimal(2) * x1.Lambda).ln()
        expected_log_factor = +(remainder.pressure.log_factor - log_divisor)

    assert axial.ordinary_correction_terms_decimal() == (expected_first, expected_second)
    assert remainder.pressure.sign != 0
    assert axial.pressure_over_two_lambda.sign == remainder.pressure.sign
    assert axial.pressure_over_two_lambda.log_scale == remainder.pressure.log_scale
    assert axial.pressure_over_two_lambda.log_factor == expected_log_factor

    # The outer Picard scale may leave the pressure far below binary64, but a
    # nonzero theorem term must still never be silently rounded to zero.
    try:
        projected = axial.pressure_over_two_lambda.to_binary64()
    except ArithmeticError as exc:
        assert "binary64" in str(exc)
    else:
        assert projected != 0.0
        assert math.isfinite(projected)


def test_zero_signed_log_pressure_stays_exactly_zero_under_picard_scaling(x1) -> None:
    # Row zero is killed by the landed regular inverse pressure chain.  The
    # outer scaling must preserve that exact zero rather than manufacture a log.
    _, axial = x1.jet_pair(0, 0, 0.0)
    assert axial.pressure_over_two_lambda.sign == 0
    assert axial.pressure_over_two_lambda.log_scale is None
    assert axial.pressure_over_two_lambda.log_factor is None


def test_factory_rejects_non_schedule_data() -> None:
    with pytest.raises(TypeError, match="data must be TailData"):
        actual_schedule_wide_first_picard_state(object(), 0.05)


def test_first_picard_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "Picard x1" in layer["capability"]
    assert "1/(2*Lambda)" in layer["capability"]
    assert "fixed-point" in layer["remaining_boundary"]
    assert "NaturalProfileAssembly" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
