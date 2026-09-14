from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
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


def _ordinary(jet):
    return jet.ordinary_factors_decimal()


def _pressure(jet):
    return jet.pressure_factors_decimal()


def test_quad1_is_bound_to_actual_x1_axisdata_and_pressure(quad1_state) -> None:
    assert isinstance(quad1_state, ActualScheduleWideFirstPicardQuad1State)
    assert quad1_state.x1.picard_x1_materialized is True
    assert quad1_state.quad1_x1_materialized is True
    assert quad1_state.natural_remainder_x1_materialized is False
    assert quad1_state.picard_x2_materialized is False
    assert quad1_state.fixed_point_materialized is False
    assert quad1_state.paper_exact is False
    assert quad1_state.Lambda == quad1_state.x1.Lambda
    assert quad1_state.epsilon == quad1_state.x1.epsilon
    assert quad1_state.angular_quadratic.epsilon == quad1_state.epsilon
    assert quad1_state.x1.remainder.axial.wide_pressure.epsilon == quad1_state.epsilon


def test_j2_row_zero_is_literal_zero_for_both_scale_families(quad1_state) -> None:
    actual = quad1_state.jet(0, 3, 0.07)
    assert isinstance(actual, MixedScaleFirstPicardQuad1CoefficientJet)
    assert _ordinary(actual) == (Decimal(0),) * 5
    assert _pressure(actual) == (Decimal(0),) * 3
    assert actual.Lambda == quad1_state.Lambda


def test_first_radial_row_matches_pinned_q_times_u1_times_phi1_over_two(quad1_state) -> None:
    eta = -0.09
    actual = quad1_state.jet(1, 0, eta)
    phi, u = quad1_state.x1.jet_pair(0, 0, eta)
    q0 = Decimal.from_float(quad1_state.angular_quadratic.jet(0, 0, eta))
    u_terms = (u.reference, u.inverse_lambda_numerator, u.inverse_lambda_squared_numerator)
    phi_terms = (phi.reference, phi.inverse_lambda_numerator, phi.inverse_lambda_squared_numerator)

    expected_ordinary = [Decimal(0) for _ in range(5)]
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for up, uv in enumerate(u_terms):
            for pp, pv in enumerate(phi_terms):
                expected_ordinary[up + pp] += q0 * uv * pv / Decimal(2)

        wide_u = +(
            quad1_state.x1.remainder.axial.pressure_normalized_factor(0, 0, eta)
            / Decimal(2)
        )
        expected_pressure = tuple(
            +(q0 * wide_u * pv / Decimal(2)) for pv in phi_terms
        )

    assert _ordinary(actual) == tuple(expected_ordinary)
    assert _pressure(actual) == expected_pressure


def test_general_jet_replays_triple_product_eta_leibniz_and_j2_exactly(quad1_state) -> None:
    n, m, eta = 3, 2, 0.04
    actual = quad1_state.jet(n, m, eta)
    input_row = n - 1
    expected_ordinary = [Decimal(0) for _ in range(5)]
    expected_pressure = [Decimal(0) for _ in range(4)]
    q = quad1_state.angular_quadratic

    with localcontext() as ctx:
        ctx.prec = PRECISION
        for qr in range(input_row + 1):
            for ur in range(input_row - qr + 1):
                pr = input_row - qr - ur
                for qm in range(m + 1):
                    remaining = m - qm
                    for um in range(remaining + 1):
                        pm = remaining - um
                        weight = Decimal(math.comb(m, qm) * math.comb(remaining, um))
                        qv = Decimal.from_float(q.jet(qr, qm, eta))
                        _, uj = quad1_state.x1.jet_pair(ur, um, eta)
                        pj, _ = quad1_state.x1.jet_pair(pr, pm, eta)
                        u_terms = (
                            uj.reference,
                            uj.inverse_lambda_numerator,
                            uj.inverse_lambda_squared_numerator,
                        )
                        p_terms = (
                            pj.reference,
                            pj.inverse_lambda_numerator,
                            pj.inverse_lambda_squared_numerator,
                        )
                        for up, uv in enumerate(u_terms):
                            for pp, pv in enumerate(p_terms):
                                expected_ordinary[up + pp] += weight * qv * uv * pv

                        wide_u = +(
                            quad1_state.x1.remainder.axial.pressure_normalized_factor(
                                ur, um, eta
                            )
                            / Decimal(2)
                        )
                        for pp, pv in enumerate(p_terms):
                            expected_pressure[1 + pp] += weight * qv * wide_u * pv

        divisor = Decimal(n) * Decimal(n + 1)
        expected_ordinary = [+(v / divisor) for v in expected_ordinary]
        expected_pressure = [+(v / divisor) for v in expected_pressure]

    assert _ordinary(actual) == tuple(expected_ordinary)
    assert _pressure(actual) == tuple(expected_pressure[1:])


def test_pressure_derived_u1_scale_survives_in_signed_log_form(quad1_state) -> None:
    candidates = [quad1_state.jet(n, 0, 0.03) for n in range(1, 5)]
    actual = next(jet for jet in candidates if any(value != 0 for value in _pressure(jet)))
    logs = actual.pressure_terms_log()
    assert len(logs) == 3
    for power, (factor, logged) in enumerate(zip(_pressure(actual), logs), start=1):
        if factor == 0:
            assert logged.sign == 0
            continue
        assert logged.sign == (1 if factor > 0 else -1)
        assert logged.log_scale == Decimal(2) * actual.amplitude_log
        with localcontext() as ctx:
            ctx.prec = PRECISION
            expected_factor = +(abs(factor).ln() - Decimal(power) * actual.Lambda.ln())
        assert logged.log_factor == expected_factor


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
    assert layer["truth_boundary"]["quad1_x1_pressure_scale_preserved"] is True
    assert layer["truth_boundary"]["natural_remainder_x1_materialized"] is False
    assert layer["truth_boundary"]["picard_x2_materialized"] is False
    assert "angularQuadraticCoefficient" in layer["capability"]
    assert "lin1/slow1" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
