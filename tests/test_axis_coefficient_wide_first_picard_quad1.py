from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_angular_square import (
    wide_first_picard_angular_square_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_quad1 import (
    ActualScheduleWideFirstPicardQuad1State,
    MixedScaleFirstPicardQuad1CoefficientJet,
    wide_first_picard_quad1_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def quad1_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_quad1_state(x1)


def _factors(jet):
    return (
        jet.reference,
        jet.inverse_lambda_numerator,
        jet.inverse_lambda_squared_numerator,
        jet.inverse_lambda_cubed_numerator,
        jet.inverse_lambda_fourth_numerator,
    )


def test_quad1_is_bound_to_genuine_first_picard_state(quad1_state) -> None:
    assert isinstance(quad1_state, ActualScheduleWideFirstPicardQuad1State)
    assert quad1_state.x1.picard_x1_materialized is True
    assert quad1_state.quad1_x1_materialized is True
    assert quad1_state.natural_remainder_x1_materialized is False
    assert quad1_state.picard_x2_materialized is False
    assert quad1_state.fixed_point_materialized is False
    assert quad1_state.paper_exact is False
    assert quad1_state.Lambda == quad1_state.x1.Lambda
    assert quad1_state.epsilon == quad1_state.x1.epsilon


def test_j2_row_zero_is_literal_zero_not_missing_coefficient_default(quad1_state) -> None:
    actual = quad1_state.jet(0, 3, 0.07)
    assert isinstance(actual, MixedScaleFirstPicardQuad1CoefficientJet)
    assert _factors(actual) == (Decimal(0),) * 5
    assert actual.Lambda == quad1_state.Lambda


def test_first_nonzero_radial_row_equals_minus_phi1_square_row_zero_over_six(
    quad1_state,
) -> None:
    eta = -0.09
    actual = quad1_state.jet(2, 0, eta)
    square = wide_first_picard_angular_square_state(quad1_state.x1).jet(0, 0, eta)

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = tuple(+(value / Decimal(-6)) for value in _factors(square))

    assert _factors(actual) == expected


def test_general_jet_replays_primitive_product_and_j2_rules_exactly(quad1_state) -> None:
    n, m, eta = 3, 2, 0.04
    actual = quad1_state.jet(n, m, eta)
    product_row = n - 1
    expected = [Decimal(0) for _ in range(5)]

    with localcontext() as ctx:
        ctx.prec = PRECISION
        for i in range(product_row + 1):
            j = product_row - i
            for k in range(m + 1):
                l = m - k
                if i == 0:
                    left = (Decimal(0), Decimal(0), Decimal(0))
                else:
                    phi_pred, _ = quad1_state.x1.jet_pair(i - 1, k, eta)
                    left = tuple(
                        +(value / Decimal(i))
                        for value in (
                            phi_pred.reference,
                            phi_pred.inverse_lambda_numerator,
                            phi_pred.inverse_lambda_squared_numerator,
                        )
                    )
                phi_right, _ = quad1_state.x1.jet_pair(j, l, eta)
                right = (
                    phi_right.reference,
                    phi_right.inverse_lambda_numerator,
                    phi_right.inverse_lambda_squared_numerator,
                )
                weight = Decimal(math.comb(m, k))
                for p, left_value in enumerate(left):
                    for q, right_value in enumerate(right):
                        expected[p + q] += weight * left_value * right_value

        divisor = Decimal(n) * Decimal(n + 1)
        expected = [+(value / -divisor) for value in expected]

    assert _factors(actual) == tuple(expected)


def test_quad1_keeps_all_inverse_lambda_families_separate(quad1_state) -> None:
    actual = quad1_state.jet(2, 1, 0.03)
    assert actual.Lambda == quad1_state.Lambda
    with localcontext() as ctx:
        ctx.prec = PRECISION
        inverse = +(Decimal(1) / actual.Lambda)
        expected = (
            +(actual.inverse_lambda_numerator * inverse),
            +(actual.inverse_lambda_squared_numerator * inverse * inverse),
            +(actual.inverse_lambda_cubed_numerator * inverse * inverse * inverse),
            +(
                actual.inverse_lambda_fourth_numerator
                * inverse
                * inverse
                * inverse
                * inverse
            ),
        )
    assert actual.correction_terms_decimal() == expected


def test_quad1_rejects_surrogate_state_and_invalid_coordinates(quad1_state) -> None:
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_quad1_state(object())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        quad1_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        quad1_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        quad1_state.jet(1, 0, 2.0)


def test_quad1_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_quad1.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert layer["truth_boundary"]["quad1_x1_materialized"] is True
    assert layer["truth_boundary"]["natural_remainder_x1_materialized"] is False
    assert layer["truth_boundary"]["picard_x2_materialized"] is False
    assert "quad1(x1)" in layer["capability"]
    assert "lin1/slow1" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
