from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_product import axis_coefficient_product
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_angular_square import (
    ActualScheduleWideFirstPicardAngularSquareState,
    MixedScaleFirstPicardAngularSquareCoefficientJet,
    wide_first_picard_angular_square_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def square_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_angular_square_state(x1)


def test_square_is_bound_to_genuine_first_picard_state(square_state) -> None:
    assert isinstance(square_state, ActualScheduleWideFirstPicardAngularSquareState)
    assert square_state.Lambda == square_state.x1.Lambda
    assert square_state.epsilon == square_state.x1.epsilon
    assert square_state.x1.picard_x1_materialized is True
    assert square_state.mixed_scale_angular_square_materialized is True
    assert square_state.natural_remainder_x1_materialized is False
    assert square_state.fixed_point_materialized is False
    assert square_state.paper_exact is False


def test_row_zero_square_has_exact_five_power_expansion(square_state) -> None:
    eta = 0.07
    phi, _ = square_state.x1.jet_pair(0, 0, eta)
    square = square_state.jet(0, 0, eta)

    assert isinstance(square, MixedScaleFirstPicardAngularSquareCoefficientJet)
    assert square.Lambda == phi.Lambda == square_state.Lambda
    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = (
            +(phi.reference * phi.reference),
            +(Decimal(2) * phi.reference * phi.inverse_lambda_numerator),
            +(
                phi.inverse_lambda_numerator * phi.inverse_lambda_numerator
                + Decimal(2)
                * phi.reference
                * phi.inverse_lambda_squared_numerator
            ),
            +(
                Decimal(2)
                * phi.inverse_lambda_numerator
                * phi.inverse_lambda_squared_numerator
            ),
            +(
                phi.inverse_lambda_squared_numerator
                * phi.inverse_lambda_squared_numerator
            ),
        )

    assert square.reference == expected[0]
    assert square.inverse_lambda_numerator == expected[1]
    assert square.inverse_lambda_squared_numerator == expected[2]
    assert square.inverse_lambda_cubed_numerator == expected[3]
    assert square.inverse_lambda_fourth_numerator == expected[4]


def test_reference_piece_cross_checks_independent_pinned_product(square_state) -> None:
    n, m, eta = 2, 1, -0.11
    square = square_state.jet(n, m, eta)
    independent = axis_coefficient_product(
        square_state.x1.reference.phi,
        square_state.x1.reference.phi,
    ).jet(n, m, eta)

    assert math.isclose(
        float(square.reference),
        independent,
        rel_tol=2e-14,
        abs_tol=2e-14,
    )


def test_all_power_numerators_follow_literal_product_convolution(square_state) -> None:
    n, m, eta = 1, 1, 0.03
    actual = square_state.jet(n, m, eta)
    expected = [Decimal(0) for _ in range(5)]

    with localcontext() as ctx:
        ctx.prec = PRECISION
        for i in range(n + 1):
            j = n - i
            for k in range(m + 1):
                l = m - k
                left, _ = square_state.x1.jet_pair(i, k, eta)
                right, _ = square_state.x1.jet_pair(j, l, eta)
                weight = Decimal(math.comb(m, k))
                for p, left_value in enumerate(
                    (
                        left.reference,
                        left.inverse_lambda_numerator,
                        left.inverse_lambda_squared_numerator,
                    )
                ):
                    for q, right_value in enumerate(
                        (
                            right.reference,
                            right.inverse_lambda_numerator,
                            right.inverse_lambda_squared_numerator,
                        )
                    ):
                        expected[p + q] += weight * left_value * right_value
        expected = [+value for value in expected]

    assert (
        actual.reference,
        actual.inverse_lambda_numerator,
        actual.inverse_lambda_squared_numerator,
        actual.inverse_lambda_cubed_numerator,
        actual.inverse_lambda_fourth_numerator,
    ) == tuple(expected)

    with localcontext() as ctx:
        ctx.prec = PRECISION
        inv = +(Decimal(1) / actual.Lambda)
        expected_corrections = (
            +(expected[1] * inv),
            +(expected[2] * inv * inv),
            +(expected[3] * inv * inv * inv),
            +(expected[4] * inv * inv * inv * inv),
        )
    assert actual.correction_terms_decimal() == expected_corrections


def test_square_rejects_invalid_indices_and_surrogate_state(square_state) -> None:
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        square_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        square_state.jet(0, True, 0.0)
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_angular_square_state(object())


def test_angular_square_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_angular_square.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "phi1 * phi1" in layer["capability"]
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert "NaturalProfileAssembly" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
