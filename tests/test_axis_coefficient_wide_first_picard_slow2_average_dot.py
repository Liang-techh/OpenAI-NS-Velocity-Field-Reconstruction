from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_data import actual_schedule_axis_coefficient_data
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_dot import (
    ActualScheduleWideFirstPicardSlow2AverageDotState,
    MixedScaleFirstPicardSlow2AverageDotCoefficientJet,
    wide_first_picard_slow2_average_dot_state,
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
    return wide_first_picard_slow2_average_dot_state(x1)


def _u_ordinary(x1, n: int, m: int, eta: float):
    _, axial = x1.jet_pair(n, m, eta)
    return (
        axial.reference,
        axial.inverse_lambda_numerator,
        axial.inverse_lambda_squared_numerator,
    )


def _u_pressure(x1, n: int, m: int, eta: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(x1.remainder.axial.pressure_normalized_factor(n, m, eta) / Decimal(2))


def _average_ordinary(x1, n: int, m: int, eta: float):
    divisor = Decimal(n + 1)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return tuple(+(value / divisor) for value in _u_ordinary(x1, n, m, eta))


def _average_pressure(x1, n: int, m: int, eta: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(_u_pressure(x1, n, m, eta) / Decimal(n + 1))


def _left_ordinary(x1, D: Decimal, n: int, m: int, eta: float):
    current = _average_ordinary(x1, n, m, eta)
    previous = _average_ordinary(x1, n, m - 1, eta) if m > 0 else None
    with localcontext() as ctx:
        ctx.prec = PRECISION
        q0 = +(Decimal(2) * D * Decimal.from_float(float(eta)))
        q1_weight = +(Decimal(m) * Decimal(2) * D)
        values = []
        for power in range(3):
            value = q0 * current[power]
            if previous is not None:
                value += q1_weight * previous[power]
            values.append(+value)
    return tuple(values)


def _left_pressure(x1, D: Decimal, n: int, m: int, eta: float) -> Decimal:
    current = _average_pressure(x1, n, m, eta)
    previous = _average_pressure(x1, n, m - 1, eta) if m > 0 else None
    with localcontext() as ctx:
        ctx.prec = PRECISION
        value = Decimal(2) * D * Decimal.from_float(float(eta)) * current
        if previous is not None:
            value += Decimal(m) * Decimal(2) * D * previous
        return +value


def _literal_dot1(branch_state, n: int, m: int, eta: float):
    if n == 0:
        return {name: Decimal(0) for name in FIELDS}

    x1 = branch_state.x1
    data = actual_schedule_axis_coefficient_data(x1.reference)
    D = Decimal.from_float(data.D)
    ordinary = [Decimal(0) for _ in range(5)]
    pressure_linear = [Decimal(0) for _ in range(4)]
    pressure_square = Decimal(0)
    source_degree = n - 1

    with localcontext() as ctx:
        ctx.prec = PRECISION
        divisor = Decimal(n * n)
        for i in range(source_degree + 1):
            j = source_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k)) * Decimal(j) / divisor
                left_o = _left_ordinary(x1, D, i, k, eta)
                right_o = _u_ordinary(x1, j, l, eta)
                left_p = _left_pressure(x1, D, i, k, eta)
                right_p = _u_pressure(x1, j, l, eta)

                for lp, lv in enumerate(left_o):
                    for rp, rv in enumerate(right_o):
                        ordinary[lp + rp] += weight * lv * rv
                for lp, lv in enumerate(left_o):
                    pressure_linear[lp + 1] += weight * lv * right_p
                for rp, rv in enumerate(right_o):
                    pressure_linear[rp + 1] += weight * left_p * rv
                pressure_square += weight * left_p * right_p

        ordinary = [+value for value in ordinary]
        pressure_linear = [+value for value in pressure_linear]
        pressure_square = +pressure_square

    return {
        "ordinary_reference": ordinary[0],
        "ordinary_inverse_lambda_numerator": ordinary[1],
        "ordinary_inverse_lambda_squared_numerator": ordinary[2],
        "ordinary_inverse_lambda_cubed_numerator": ordinary[3],
        "ordinary_inverse_lambda_fourth_numerator": ordinary[4],
        "pressure_linear_inverse_lambda_numerator": pressure_linear[1],
        "pressure_linear_inverse_lambda_squared_numerator": pressure_linear[2],
        "pressure_linear_inverse_lambda_cubed_numerator": pressure_linear[3],
        "pressure_square_inverse_lambda_squared_numerator": pressure_square,
    }


def _literal_reference_only(branch_state, n: int, m: int, eta: float) -> Decimal:
    if n == 0:
        return Decimal(0)
    x1 = branch_state.x1
    data = actual_schedule_axis_coefficient_data(x1.reference)
    D = Decimal.from_float(data.D)
    source_degree = n - 1
    total = Decimal(0)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        divisor = Decimal(n * n)
        eta_decimal = Decimal.from_float(float(eta))
        for i in range(source_degree + 1):
            j = source_degree - i
            for k in range(m + 1):
                l = m - k
                current = Decimal.from_float(x1.reference.u.jet(i, k, eta)) / Decimal(i + 1)
                previous = (
                    Decimal.from_float(x1.reference.u.jet(i, k - 1, eta)) / Decimal(i + 1)
                    if k > 0
                    else None
                )
                left = Decimal(2) * D * eta_decimal * current
                if previous is not None:
                    left += Decimal(k) * Decimal(2) * D * previous
                right = Decimal.from_float(x1.reference.u.jet(j, l, eta))
                total += Decimal(math.comb(m, k)) * Decimal(j) / divisor * left * right
        return +total


def test_branch_is_bound_to_genuine_first_picard_state(branch_state) -> None:
    assert isinstance(branch_state, ActualScheduleWideFirstPicardSlow2AverageDotState)
    assert branch_state.Lambda == branch_state.x1.Lambda
    assert branch_state.epsilon == branch_state.x1.epsilon
    assert branch_state.x1.picard_x1_materialized is True
    assert branch_state.slow2_average_dot_branch_materialized is True
    assert branch_state.slow2_complete is False
    assert branch_state.natural_remainder_x1_materialized is False
    assert branch_state.fixed_point_materialized is False
    assert branch_state.paper_exact is False


def test_ordinary_reference_follows_independent_reference_only_dot1_formula(branch_state) -> None:
    n, m, eta = 4, 2, -0.11
    actual = branch_state.jet(n, m, eta)
    expected = _literal_reference_only(branch_state, n, m, eta)
    assert actual.ordinary_reference == expected


def test_all_scale_numerators_follow_literal_dot1_formula(branch_state) -> None:
    n, m, eta = 4, 2, 0.03
    actual = branch_state.jet(n, m, eta)
    literal = _literal_dot1(branch_state, n, m, eta)

    assert isinstance(actual, MixedScaleFirstPicardSlow2AverageDotCoefficientJet)
    for name in FIELDS:
        assert getattr(actual, name) == literal[name]

    assert actual.Lambda == branch_state.Lambda
    assert actual.amplitude_log == (
        branch_state.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
    )

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


def test_dot1_row_zero_is_exact_zero_without_skipping_upstream_state(branch_state) -> None:
    actual = branch_state.jet(0, 2, 0.07)
    for name in FIELDS:
        assert getattr(actual, name) == 0
    assert actual.Lambda == branch_state.Lambda


def test_nonzero_pressure_terms_keep_a2_and_a4_scales(branch_state) -> None:
    eta = 0.07
    selected = next(
        (
            jet
            for n in range(2, 13)
            for m in range(0, 3)
            for jet in (branch_state.jet(n, m, eta),)
            if any(
                getattr(jet, name) != 0
                for name in (
                    "pressure_linear_inverse_lambda_numerator",
                    "pressure_linear_inverse_lambda_squared_numerator",
                    "pressure_linear_inverse_lambda_cubed_numerator",
                )
            )
            and jet.pressure_square_inverse_lambda_squared_numerator != 0
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
        for power, (numerator, logged) in enumerate(
            zip(linear_numerators, actual.pressure_linear_terms_log()),
            start=1,
        ):
            if numerator == 0:
                assert logged.sign == 0
                continue
            assert logged.sign == (1 if numerator > 0 else -1)
            assert logged.log_scale == expected_a2
            assert logged.log_factor == +(
                abs(numerator).ln() - Decimal(power) * lambda_log
            )

        square_numerator = actual.pressure_square_inverse_lambda_squared_numerator
        square_log = actual.pressure_square_term_log()
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
        wide_first_picard_slow2_average_dot_state(object())


def test_branch_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_average_dot.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "dot1" in layer["capability"]
    assert "2*D*eta" in layer["capability"]
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert "NaturalProfileAssembly" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
