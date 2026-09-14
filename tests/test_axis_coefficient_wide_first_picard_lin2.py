from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_lin2 import (
    ActualScheduleWideFirstPicardLin2State,
    MixedScaleFirstPicardLin2CoefficientJet,
    wide_first_picard_lin2_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def lin2_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_lin2_state(x1)


def _decimal(value: float) -> Decimal:
    value = float(value)
    assert math.isfinite(value)
    return Decimal.from_float(value)


def _product_jet(left, right, n: int, m: int, eta: float) -> Decimal:
    """Independent replay of the landed binary64 coefficient-product semantics."""

    terms = []
    for i in range(n + 1):
        j = n - i
        for k in range(m + 1):
            l = m - k
            terms.append(
                float(math.comb(m, k))
                * left.jet(i, k, eta)
                * right.jet(j, l, eta)
            )
    return _decimal(math.fsum(terms))


def _axial_linear_jet(state, n: int, m: int, eta: float) -> Decimal:
    d = state.axis_data
    with localcontext() as ctx:
        ctx.prec = PRECISION
        A = _decimal(d.A)
        return +(
            A * _decimal(d.one.jet(n, m, eta))
            - Decimal(4) * A * _product_jet(d.eta, d.uStar, n, m, eta)
            + _product_jet(d.d, d.uStarEta, n, m, eta)
        )


def _u_ordinary_terms(state, n: int, m: int, eta: float):
    _, u = state.x1.jet_pair(n, m, eta)
    return (
        u.reference,
        u.inverse_lambda_numerator,
        u.inverse_lambda_squared_numerator,
    )


def _u_pressure_term(state, n: int, m: int, eta: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(
            state.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
            / Decimal(2)
        )


def test_lin2_is_bound_to_genuine_x1_and_actual_axis_data(lin2_state) -> None:
    assert isinstance(lin2_state, ActualScheduleWideFirstPicardLin2State)
    assert lin2_state.x1.picard_x1_materialized is True
    assert lin2_state.lin2_x1_materialized is True
    assert lin2_state.natural_remainder_x1_materialized is False
    assert lin2_state.picard_x2_materialized is False
    assert lin2_state.fixed_point_materialized is False
    assert lin2_state.paper_exact is False
    assert lin2_state.Lambda == lin2_state.x1.Lambda
    assert lin2_state.epsilon == lin2_state.x1.epsilon
    assert lin2_state.axis_data.epsilon == lin2_state.epsilon
    assert lin2_state.operators.epsilon == lin2_state.epsilon


def test_all_genuine_lin2_scales_have_literal_zero_row(lin2_state) -> None:
    actual = lin2_state.jet(0, 4, 0.07)
    assert isinstance(actual, MixedScaleFirstPicardLin2CoefficientJet)
    assert actual.ordinary_factors_decimal() == (
        Decimal(0),
        Decimal(0),
        Decimal(0),
    )
    assert actual.pressure_factor_decimal() == Decimal(0)
    assert actual.Lambda == lin2_state.Lambda


def test_first_radial_row_matches_j1_plus_param1_identity(lin2_state) -> None:
    n, m, eta = 1, 0, -0.09
    actual = lin2_state.jet(n, m, eta)
    fixed = _axial_linear_jet(lin2_state, 0, 0, eta)
    hstar = _decimal(lin2_state.axis_data.hStar.jet(0, 0, eta))

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected_ordinary = tuple(
            +(fixed * value + derivative * hstar)
            for value, derivative in zip(
                _u_ordinary_terms(lin2_state, 0, 0, eta),
                _u_ordinary_terms(lin2_state, 0, 1, eta),
            )
        )
        expected_pressure = +(
            fixed * _u_pressure_term(lin2_state, 0, 0, eta)
            + _u_pressure_term(lin2_state, 0, 1, eta) * hstar
        )

    # At output row one, dot1 has second-argument radial degree zero.
    assert actual.ordinary_factors_decimal() == expected_ordinary
    assert actual.pressure_factor_decimal() == expected_pressure


def test_general_jet_replays_all_three_pinned_branches_exactly(lin2_state) -> None:
    n, m, eta = 3, 2, 0.04
    actual = lin2_state.jet(n, m, eta)
    d = lin2_state.axis_data
    radial_degree = n - 1
    ordinary = [Decimal(0), Decimal(0), Decimal(0)]
    pressure = Decimal(0)

    with localcontext() as ctx:
        ctx.prec = PRECISION

        # j1(product(axialLinearCoefficient(d), u1))
        for i in range(radial_degree + 1):
            j = radial_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k))
                fixed = _axial_linear_jet(lin2_state, i, k, eta)
                for power, u_value in enumerate(
                    _u_ordinary_terms(lin2_state, j, l, eta)
                ):
                    ordinary[power] += weight * fixed * u_value
                pressure += (
                    weight
                    * fixed
                    * _u_pressure_term(lin2_state, j, l, eta)
                )

        # dot1(wStar, u1)
        for i in range(radial_degree + 1):
            j = radial_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k) * j)
                fixed = _decimal(d.wStar.jet(i, k, eta))
                for power, u_value in enumerate(
                    _u_ordinary_terms(lin2_state, j, l, eta)
                ):
                    ordinary[power] += weight * fixed * u_value
                pressure += (
                    weight
                    * fixed
                    * _u_pressure_term(lin2_state, j, l, eta)
                )

        # param1(u1, hStar)
        for i in range(radial_degree + 1):
            j = radial_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k))
                fixed = _decimal(d.hStar.jet(j, l, eta))
                for power, u_value in enumerate(
                    _u_ordinary_terms(lin2_state, i, k + 1, eta)
                ):
                    ordinary[power] += weight * u_value * fixed
                pressure += (
                    weight
                    * _u_pressure_term(lin2_state, i, k + 1, eta)
                    * fixed
                )

        divisor = Decimal(n) * Decimal(n)
        expected_ordinary = tuple(+(value / divisor) for value in ordinary)
        expected_pressure = +(pressure / divisor)

    assert actual.ordinary_factors_decimal() == expected_ordinary
    assert actual.pressure_factor_decimal() == expected_pressure


def test_lin2_preserves_genuine_pressure_family_without_binary64_collapse(
    lin2_state,
) -> None:
    actual = lin2_state.jet(2, 1, 0.03)
    assert len(actual.ordinary_factors_decimal()) == 3
    assert isinstance(actual.pressure_factor_decimal(), Decimal)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected_log = +(
            lin2_state.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(0.03)
        )
    assert actual.amplitude_log == expected_log
    signed_log = actual.pressure_term_log()
    assert signed_log.sign in (-1, 0, 1)


def test_lin2_rejects_surrogate_state_and_invalid_coordinates(lin2_state) -> None:
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_lin2_state(object())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        lin2_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        lin2_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        lin2_state.jet(1, 0, 2.0)


def test_lin2_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_lin2.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert layer["truth_boundary"]["lin1_x1_materialized"] is True
    assert layer["truth_boundary"]["slow1_x1_materialized"] is False
    assert layer["truth_boundary"]["lin2_x1_materialized"] is True
    assert layer["truth_boundary"]["natural_remainder_x1_materialized"] is False
    assert layer["truth_boundary"]["picard_x2_materialized"] is False
    assert "axialLinearCoefficient" in layer["capability"]
    assert "a^2/Lambda" in layer["capability"]
    assert "slow1(x1)" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
