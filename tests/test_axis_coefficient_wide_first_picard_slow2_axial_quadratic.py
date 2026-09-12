from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_data import actual_schedule_axis_coefficient_data
from openai_ns_reconstruction.axis_coefficient_product import axis_coefficient_product
from openai_ns_reconstruction.axis_coefficient_regular_inverse import axis_coefficient_regular_inverse
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_axial_square import (
    wide_first_picard_axial_square_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_axial_quadratic import (
    ActualScheduleWideFirstPicardSlow2AxialQuadraticState,
    MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet,
    wide_first_picard_slow2_axial_quadratic_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96
FIELDS = (
    "ordinary_reference",
    "ordinary_inverse_lambda_numerator",
    "ordinary_inverse_lambda_squared_numerator",
    "ordinary_inverse_lambda_cubed_numerator",
    "ordinary_inverse_lambda_fourth_numerator",
    "pressure_linear_inverse_lambda_numerator",
    "pressure_linear_inverse_lambda_squared_numerator",
    "pressure_linear_inverse_lambda_cubed_numerator",
    "pressure_square_inverse_lambda_squared_numerator",
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def branch_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_slow2_axial_quadratic_state(x1)


def _literal_branch(branch_state, n: int, m: int, eta: float):
    square = wide_first_picard_axial_square_state(branch_state.x1)
    if n == 0:
        source = square.jet(0, m, eta)
        return {name: Decimal(0) for name in FIELDS}, source

    current = square.jet(n - 1, m, eta)
    previous = square.jet(n - 1, m - 1, eta) if m > 0 else None
    data = actual_schedule_axis_coefficient_data(branch_state.x1.reference)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        A = Decimal.from_float(data.A)
        q0 = +(Decimal(2) * A * Decimal.from_float(float(eta)))
        q1_weight = +(Decimal(m) * Decimal(2) * A)
        divisor = Decimal(n * n)
        values = {}
        for name in FIELDS:
            value = q0 * getattr(current, name)
            if previous is not None:
                value += q1_weight * getattr(previous, name)
            values[name] = +(value / divisor)
    return values, current


def test_branch_is_bound_to_genuine_first_picard_state(branch_state) -> None:
    assert isinstance(branch_state, ActualScheduleWideFirstPicardSlow2AxialQuadraticState)
    assert branch_state.Lambda == branch_state.x1.Lambda
    assert branch_state.epsilon == branch_state.x1.epsilon
    assert branch_state.x1.picard_x1_materialized is True
    assert branch_state.slow2_axial_quadratic_branch_materialized is True
    assert branch_state.slow2_complete is False
    assert branch_state.natural_remainder_x1_materialized is False
    assert branch_state.fixed_point_materialized is False
    assert branch_state.paper_exact is False


def test_ordinary_reference_cross_checks_independent_pinned_operator_path(branch_state) -> None:
    n, m, eta = 3, 1, -0.11
    actual = branch_state.jet(n, m, eta)
    data = actual_schedule_axis_coefficient_data(branch_state.x1.reference)
    reference_square = axis_coefficient_product(
        branch_state.x1.reference.u,
        branch_state.x1.reference.u,
    )
    eta_times_square = axis_coefficient_product(data.eta, reference_square)
    independent = axis_coefficient_regular_inverse(eta_times_square, 1).jet(n, m, eta)
    independent *= 2.0 * data.A

    assert math.isclose(
        float(actual.ordinary_reference),
        independent,
        rel_tol=3e-13,
        abs_tol=3e-13,
    )


def test_all_scale_numerators_follow_literal_eta_product_and_j1(branch_state) -> None:
    n, m, eta = 3, 2, 0.03
    actual = branch_state.jet(n, m, eta)
    literal, source = _literal_branch(branch_state, n, m, eta)

    assert isinstance(actual, MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet)
    for name in FIELDS:
        assert getattr(actual, name) == literal[name]
    assert actual.Lambda == source.Lambda
    assert actual.amplitude_log == source.amplitude_log

    with localcontext() as ctx:
        ctx.prec = PRECISION
        inverse = +(Decimal(1) / actual.Lambda)
        expected = (
            +(literal["ordinary_inverse_lambda_numerator"] * inverse),
            +(literal["ordinary_inverse_lambda_squared_numerator"] * inverse * inverse),
            +(literal["ordinary_inverse_lambda_cubed_numerator"] * inverse * inverse * inverse),
            +(literal["ordinary_inverse_lambda_fourth_numerator"] * inverse * inverse * inverse * inverse),
        )
    assert actual.ordinary_correction_terms_decimal() == expected


def test_j1_row_zero_is_exact_zero_without_skipping_upstream_state(branch_state) -> None:
    actual = branch_state.jet(0, 2, 0.07)
    for name in FIELDS:
        assert getattr(actual, name) == 0
    assert actual.Lambda == branch_state.Lambda


def test_nonzero_pressure_terms_keep_a2_and_a4_scales(branch_state) -> None:
    eta = 0.07
    selected = next(
        (
            branch_state.jet(n, m, eta)
            for n in range(1, 12)
            for m in range(0, 3)
            if any(
                getattr(branch_state.jet(n, m, eta), name) != 0
                for name in (
                    "pressure_linear_inverse_lambda_numerator",
                    "pressure_linear_inverse_lambda_squared_numerator",
                    "pressure_linear_inverse_lambda_cubed_numerator",
                )
            )
            and branch_state.jet(n, m, eta).pressure_square_inverse_lambda_squared_numerator != 0
        ),
        None,
    )
    assert selected is not None
    actual = selected

    with localcontext() as ctx:
        ctx.prec = PRECISION
        lambda_log = +actual.Lambda.ln()
        expected_a2 = +(Decimal(2) * actual.amplitude_log)
        expected_a4 = +(Decimal(4) * actual.amplitude_log)
        linear_numerators = (
            actual.pressure_linear_inverse_lambda_numerator,
            actual.pressure_linear_inverse_lambda_squared_numerator,
            actual.pressure_linear_inverse_lambda_cubed_numerator,
        )
        linear_logs = actual.pressure_linear_terms_log()
        for power, (numerator, logged) in enumerate(zip(linear_numerators, linear_logs), start=1):
            if numerator == 0:
                assert logged.sign == 0
                continue
            assert logged.sign == (1 if numerator > 0 else -1)
            assert logged.log_scale == expected_a2
            assert logged.log_factor == +(
                abs(numerator).ln() - Decimal(power) * lambda_log
            )

        square_log = actual.pressure_square_term_log()
        square_numerator = actual.pressure_square_inverse_lambda_squared_numerator
        assert square_log.sign == (1 if square_numerator > 0 else -1)
        assert square_log.log_scale == expected_a4
        assert square_log.log_factor == +(
            abs(square_numerator).ln() - Decimal(2) * lambda_log
        )

    nonzero_linear = next(term for term in actual.pressure_linear_terms_log() if term.sign != 0)
    with pytest.raises(ArithmeticError, match="binary64"):
        nonzero_linear.to_binary64()
    with pytest.raises(ArithmeticError, match="binary64"):
        actual.pressure_square_term_log().to_binary64()


def test_branch_rejects_invalid_indices_window_and_surrogate_state(branch_state) -> None:
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        branch_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        branch_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        branch_state.jet(1, 0, 2.0)
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_slow2_axial_quadratic_state(object())


def test_branch_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_axial_quadratic.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "j1" in layer["capability"]
    assert "2*A*eta" in layer["capability"]
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert "NaturalProfileAssembly" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
