from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_lin1 import (
    ActualScheduleWideFirstPicardLin1State,
    MixedScaleFirstPicardLin1CoefficientJet,
    wide_first_picard_lin1_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def lin1_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_lin1_state(x1)


def _factors(jet):
    return jet.factors_decimal()


def _decimal(value: float) -> Decimal:
    value = float(value)
    assert math.isfinite(value)
    return Decimal.from_float(value)


def _product_jet(left, right, n: int, m: int, eta: float) -> Decimal:
    """Independent exact-Decimal replay of the landed coefficient product."""

    out = Decimal(0)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for i in range(n + 1):
            j = n - i
            for k in range(m + 1):
                l = m - k
                out += (
                    Decimal(math.comb(m, k))
                    * _decimal(left.jet(i, k, eta))
                    * _decimal(right.jet(j, l, eta))
                )
        return +out


def _angular_linear_jet(state, n: int, m: int, eta: float) -> Decimal:
    d = state.axis_data
    with localcontext() as ctx:
        ctx.prec = PRECISION
        h = _decimal(d.h)
        return +(
            _decimal(d.wStar.jet(n, m, eta))
            + h * _decimal(d.one.jet(n, m, eta))
            - Decimal(2) * h * _product_jet(d.eta, d.uStar, n, m, eta)
        )


def _phi_terms(state, n: int, m: int, eta: float):
    phi, _ = state.x1.jet_pair(n, m, eta)
    return (
        phi.reference,
        phi.inverse_lambda_numerator,
        phi.inverse_lambda_squared_numerator,
    )


def test_lin1_is_bound_to_genuine_x1_and_actual_axis_data(lin1_state) -> None:
    assert isinstance(lin1_state, ActualScheduleWideFirstPicardLin1State)
    assert lin1_state.x1.picard_x1_materialized is True
    assert lin1_state.lin1_x1_materialized is True
    assert lin1_state.natural_remainder_x1_materialized is False
    assert lin1_state.picard_x2_materialized is False
    assert lin1_state.fixed_point_materialized is False
    assert lin1_state.paper_exact is False
    assert lin1_state.Lambda == lin1_state.x1.Lambda
    assert lin1_state.epsilon == lin1_state.x1.epsilon
    assert lin1_state.axis_data.epsilon == lin1_state.epsilon
    assert lin1_state.operators.epsilon == lin1_state.epsilon


def test_all_three_regular_inverse_branches_have_literal_zero_row(lin1_state) -> None:
    actual = lin1_state.jet(0, 4, 0.07)
    assert isinstance(actual, MixedScaleFirstPicardLin1CoefficientJet)
    assert _factors(actual) == (Decimal(0), Decimal(0), Decimal(0))
    assert actual.Lambda == lin1_state.Lambda


def test_first_radial_row_matches_j2_plus_param2_identity(lin1_state) -> None:
    n, m, eta = 1, 0, -0.09
    actual = lin1_state.jet(n, m, eta)
    angular_linear = _angular_linear_jet(lin1_state, 0, 0, eta)
    hstar = _decimal(lin1_state.axis_data.hStar.jet(0, 0, eta))
    phi_value = _phi_terms(lin1_state, 0, 0, eta)
    phi_eta = _phi_terms(lin1_state, 0, 1, eta)

    # At output row one, the dot2 Euler factor is j=0 and vanishes exactly.
    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = tuple(
            +(angular_linear * value + derivative * hstar) / Decimal(2)
            for value, derivative in zip(phi_value, phi_eta)
        )
    assert _factors(actual) == expected


def test_general_jet_replays_all_three_pinned_branches_exactly(lin1_state) -> None:
    n, m, eta = 3, 2, 0.04
    actual = lin1_state.jet(n, m, eta)
    d = lin1_state.axis_data
    radial_degree = n - 1
    expected = [Decimal(0), Decimal(0), Decimal(0)]

    with localcontext() as ctx:
        ctx.prec = PRECISION
        # j2(product(angularLinearCoefficient(d), phi1))
        for i in range(radial_degree + 1):
            j = radial_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k))
                fixed = _angular_linear_jet(lin1_state, i, k, eta)
                for power, phi_value in enumerate(_phi_terms(lin1_state, j, l, eta)):
                    expected[power] += weight * fixed * phi_value

        # dot2(wStar, phi1)
        for i in range(radial_degree + 1):
            j = radial_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k) * j)
                fixed = _decimal(d.wStar.jet(i, k, eta))
                for power, phi_value in enumerate(_phi_terms(lin1_state, j, l, eta)):
                    expected[power] += weight * fixed * phi_value

        # param2(phi1, hStar)
        for i in range(radial_degree + 1):
            j = radial_degree - i
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k))
                fixed = _decimal(d.hStar.jet(j, l, eta))
                for power, phi_value in enumerate(_phi_terms(lin1_state, i, k + 1, eta)):
                    expected[power] += weight * phi_value * fixed

        divisor = Decimal(n) * Decimal(n + 1)
        expected = tuple(+(value / divisor) for value in expected)

    assert _factors(actual) == expected


def test_lin1_preserves_only_the_three_genuine_phi1_scales(lin1_state) -> None:
    actual = lin1_state.jet(2, 1, 0.03)
    assert len(_factors(actual)) == 3
    assert not hasattr(actual, "pressure_inverse_lambda_numerator")
    corrections = actual.correction_terms_decimal()
    assert len(corrections) == 2
    with localcontext() as ctx:
        ctx.prec = PRECISION
        assert corrections[0] == +(actual.inverse_lambda_numerator / actual.Lambda)
        assert corrections[1] == +(
            actual.inverse_lambda_squared_numerator / actual.Lambda / actual.Lambda
        )


def test_lin1_rejects_surrogate_state_and_invalid_coordinates(lin1_state) -> None:
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_lin1_state(object())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        lin1_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        lin1_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        lin1_state.jet(1, 0, 2.0)


def test_lin1_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_lin1.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert layer["truth_boundary"]["lin1_x1_materialized"] is True
    assert layer["truth_boundary"]["slow1_x1_materialized"] is False
    assert layer["truth_boundary"]["lin2_x1_materialized"] is False
    assert layer["truth_boundary"]["natural_remainder_x1_materialized"] is False
    assert layer["truth_boundary"]["picard_x2_materialized"] is False
    assert "angularLinearCoefficient" in layer["capability"]
    assert "slow1(x1)" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
