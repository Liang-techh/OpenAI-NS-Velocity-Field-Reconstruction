from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_product import axis_coefficient_product
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_axial_square import (
    ActualScheduleWideFirstPicardAxialSquareState,
    MixedScaleFirstPicardAxialSquareCoefficientJet,
    wide_first_picard_axial_square_state,
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
    return wide_first_picard_axial_square_state(x1)


def _pressure_numerator(square_state, n: int, m: int, eta: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(
            square_state.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
            / Decimal(2)
        )


def _literal_square(square_state, n: int, m: int, eta: float):
    ordinary = [Decimal(0) for _ in range(5)]
    pressure_linear = [Decimal(0) for _ in range(4)]
    pressure_square = Decimal(0)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for i in range(n + 1):
            j = n - i
            for k in range(m + 1):
                l = m - k
                _, left = square_state.x1.jet_pair(i, k, eta)
                _, right = square_state.x1.jet_pair(j, l, eta)
                weight = Decimal(math.comb(m, k))
                left_ordinary = (
                    left.reference,
                    left.inverse_lambda_numerator,
                    left.inverse_lambda_squared_numerator,
                )
                right_ordinary = (
                    right.reference,
                    right.inverse_lambda_numerator,
                    right.inverse_lambda_squared_numerator,
                )
                for p, left_value in enumerate(left_ordinary):
                    for q, right_value in enumerate(right_ordinary):
                        ordinary[p + q] += weight * left_value * right_value

                left_pressure = _pressure_numerator(square_state, i, k, eta)
                right_pressure = _pressure_numerator(square_state, j, l, eta)
                for p, left_value in enumerate(left_ordinary):
                    pressure_linear[p + 1] += weight * left_value * right_pressure
                for q, right_value in enumerate(right_ordinary):
                    pressure_linear[q + 1] += weight * left_pressure * right_value
                pressure_square += weight * left_pressure * right_pressure

        return (
            tuple(+value for value in ordinary),
            tuple(+value for value in pressure_linear),
            +pressure_square,
        )


def test_axial_square_is_bound_to_genuine_first_picard_state(square_state) -> None:
    assert isinstance(square_state, ActualScheduleWideFirstPicardAxialSquareState)
    assert square_state.Lambda == square_state.x1.Lambda
    assert square_state.epsilon == square_state.x1.epsilon
    assert square_state.x1.picard_x1_materialized is True
    assert square_state.mixed_scale_axial_square_materialized is True
    assert square_state.signed_log_pressure_product_materialized is True
    assert square_state.natural_remainder_x1_materialized is False
    assert square_state.fixed_point_materialized is False
    assert square_state.paper_exact is False


def test_ordinary_reference_piece_cross_checks_independent_pinned_product(square_state) -> None:
    n, m, eta = 2, 1, -0.11
    square = square_state.jet(n, m, eta)
    independent = axis_coefficient_product(
        square_state.x1.reference.u,
        square_state.x1.reference.u,
    ).jet(n, m, eta)

    assert math.isclose(
        float(square.ordinary_reference),
        independent,
        rel_tol=2e-14,
        abs_tol=2e-14,
    )


def test_all_scale_numerators_follow_literal_product_convolution(square_state) -> None:
    n, m, eta = 2, 1, 0.03
    actual = square_state.jet(n, m, eta)
    ordinary, pressure_linear, pressure_square = _literal_square(
        square_state, n, m, eta
    )

    assert isinstance(actual, MixedScaleFirstPicardAxialSquareCoefficientJet)
    assert (
        actual.ordinary_reference,
        actual.ordinary_inverse_lambda_numerator,
        actual.ordinary_inverse_lambda_squared_numerator,
        actual.ordinary_inverse_lambda_cubed_numerator,
        actual.ordinary_inverse_lambda_fourth_numerator,
    ) == ordinary
    assert (
        actual.pressure_linear_inverse_lambda_numerator,
        actual.pressure_linear_inverse_lambda_squared_numerator,
        actual.pressure_linear_inverse_lambda_cubed_numerator,
    ) == pressure_linear[1:]
    assert actual.pressure_square_inverse_lambda_squared_numerator == pressure_square

    with localcontext() as ctx:
        ctx.prec = PRECISION
        inverse = +(Decimal(1) / actual.Lambda)
        expected_ordinary_corrections = (
            +(ordinary[1] * inverse),
            +(ordinary[2] * inverse * inverse),
            +(ordinary[3] * inverse * inverse * inverse),
            +(ordinary[4] * inverse * inverse * inverse * inverse),
        )
    assert actual.ordinary_correction_terms_decimal() == expected_ordinary_corrections


def test_nonzero_pressure_terms_keep_a2_and_a4_log_scales(square_state) -> None:
    eta = 0.0
    selected = next(
        (
            (row, literal)
            for row in range(1, 10)
            for literal in [_literal_square(square_state, row, 0, eta)]
            if any(value != 0 for value in literal[1][1:])
            and literal[2] != 0
        ),
        None,
    )
    assert selected is not None
    n, (ordinary, pressure_linear, pressure_square) = selected

    actual = square_state.jet(n, 0, eta)
    assert any(value != 0 for value in pressure_linear[1:])
    assert pressure_square != 0
    assert ordinary[0] == actual.ordinary_reference

    linear_logs = actual.pressure_linear_terms_log()
    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected_a2_scale = +(Decimal(2) * actual.amplitude_log)
        expected_a4_scale = +(Decimal(4) * actual.amplitude_log)
        lambda_log = +actual.Lambda.ln()
        for power, (numerator, logged) in enumerate(
            zip(pressure_linear[1:], linear_logs), start=1
        ):
            if numerator == 0:
                assert logged.sign == 0
                continue
            assert logged.sign == (1 if numerator > 0 else -1)
            assert logged.log_scale == expected_a2_scale
            assert logged.log_factor == +(
                abs(numerator).ln() - Decimal(power) * lambda_log
            )

        square_log = actual.pressure_square_term_log()
        assert square_log.sign == (1 if pressure_square > 0 else -1)
        assert square_log.log_scale == expected_a4_scale
        assert square_log.log_factor == +(
            abs(pressure_square).ln() - Decimal(2) * lambda_log
        )

    # The theorem-selected amplitude is currently far below binary64 range.
    # A nonzero product contribution must therefore fail closed, not become 0.
    nonzero_linear = next(term for term in linear_logs if term.sign != 0)
    with pytest.raises(ArithmeticError, match="binary64"):
        nonzero_linear.to_binary64()
    with pytest.raises(ArithmeticError, match="binary64"):
        actual.pressure_square_term_log().to_binary64()


def test_axial_square_rejects_invalid_indices_and_surrogate_state(square_state) -> None:
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        square_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        square_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        square_state.jet(0, 0, 2.0)
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_axial_square_state(object())


def test_axial_square_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_axial_square.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "u1 * u1" in layer["capability"]
    assert "a^2" in layer["capability"]
    assert "a^4" in layer["capability"]
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert "NaturalProfileAssembly" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
