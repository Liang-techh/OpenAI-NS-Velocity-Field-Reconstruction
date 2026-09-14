from decimal import Decimal, localcontext
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow1 import (
    ActualScheduleWideFirstPicardSlow1State,
    wide_first_picard_slow1_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData

PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def state():
    return wide_first_picard_slow1_state(
        actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    )


def _conv(a, b):
    out = [Decimal(0) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def _phi(state, n, m, eta):
    p, _ = state.x1.jet_pair(n, m, eta)
    return [p.reference, p.inverse_lambda_numerator, p.inverse_lambda_squared_numerator]


def _u(state, n, m, eta):
    _, u = state.x1.jet_pair(n, m, eta)
    return [u.reference, u.inverse_lambda_numerator, u.inverse_lambda_squared_numerator]


def _up(state, n, m, eta):
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(
            state.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
            / Decimal(2)
        )


def _bu(state, n, m, eta):
    return [value / Decimal(n + 1) for value in _u(state, n, m, eta)]


def _bup(state, n, m, eta):
    return _up(state, n, m, eta) / Decimal(n + 1)


def test_slow1_is_bound_to_genuine_x1_and_remains_fail_closed(state):
    assert isinstance(state, ActualScheduleWideFirstPicardSlow1State)
    assert state.x1.picard_x1_materialized is True
    assert state.slow1_x1_materialized is True
    assert state.natural_remainder_x1_materialized is False
    assert state.picard_x2_materialized is False
    assert state.fixed_point_materialized is False
    assert state.paper_exact is False
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_slow1_state(object())


def test_pinned_average_is_exact_row_divisor(state):
    n, m, eta = 3, 2, 0.03
    assert state._bu(n, m, eta) == tuple(
        value / Decimal(n + 1) for value in _u(state, n, m, eta)
    )
    assert state._bup(n, m, eta) == _up(state, n, m, eta) / Decimal(n + 1)


def test_slow1_has_literal_zero_inverse_row(state):
    actual = state.jet(0, 3, -0.02)
    assert actual.ordinary_factors_decimal() == (Decimal(0),) * 5
    assert actual.pressure_factors_decimal() == (Decimal(0),) * 3


def test_first_radial_row_replays_pinned_identity_exactly(state):
    # At n=1 both dot2 and mixed2 vanish because their second radial degree is 0.
    eta = -0.07
    actual = state.jet(1, 0, eta)
    d = state.axis_data
    D = Decimal.from_float(float(d.D))
    h = Decimal.from_float(float(d.h))
    eta0 = Decimal.from_float(float(d.eta.jet(0, 0, eta)))
    d0 = Decimal.from_float(float(d.d.jet(0, 0, eta)))

    with localcontext() as ctx:
        ctx.prec = PRECISION
        avgc = Decimal(2) * D * eta0
        slowc = Decimal(2) * h * eta0
        left_o = [
            avgc * a + slowc * b
            for a, b in zip(_bu(state, 0, 0, eta), _u(state, 0, 0, eta))
        ]
        left_p = avgc * _bup(state, 0, 0, eta) + slowc * _up(state, 0, 0, eta)
        phi0 = _phi(state, 0, 0, eta)
        ordinary = _conv(left_o, phi0)
        pressure = [Decimal(0)] * 4
        for q, value in enumerate(phi0):
            pressure[q + 1] += left_p * value

        dphi0 = [d0 * value for value in phi0]
        ordinary = [
            a + b for a, b in zip(ordinary, _conv(_bu(state, 0, 1, eta), dphi0))
        ]
        for q, value in enumerate(dphi0):
            pressure[q + 1] += _bup(state, 0, 1, eta) * value

        du0 = [d0 * value for value in _u(state, 0, 0, eta)]
        ordinary = [
            a - b for a, b in zip(ordinary, _conv(_phi(state, 0, 1, eta), du0))
        ]
        for q, value in enumerate(_phi(state, 0, 1, eta)):
            pressure[q + 1] -= value * d0 * _up(state, 0, 0, eta)

        expected_o = tuple(+(value / Decimal(2)) for value in ordinary)
        expected_p = tuple(+(pressure[q] / Decimal(2)) for q in range(1, 4))

    assert actual.ordinary_factors_decimal() == expected_o
    assert actual.pressure_factors_decimal() == expected_p


def test_slow1_preserves_pressure_linear_families_and_validates_coordinates(state):
    actual = state.jet(3, 2, 0.04)
    assert len(actual.ordinary_factors_decimal()) == 5
    assert len(actual.pressure_factors_decimal()) == 3
    assert all(isinstance(value, Decimal) for value in actual.pressure_factors_decimal())
    assert all(term.sign in (-1, 0, 1) for term in actual.pressure_terms_log())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        state.jet(1, 0, 2.0)


def test_slow1_provenance_is_machine_readable_and_fail_closed():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow1.json"
        ).read_text(encoding="utf-8")
    )
    layer = manifest["layer"]
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    assert layer["status"] == "formal-structure"
    assert layer["truth_boundary"]["slow1_x1_materialized"] is True
    assert layer["truth_boundary"]["natural_remainder_x1_materialized"] is False
    assert layer["truth_boundary"]["picard_x2_materialized"] is False
    assert layer["truth_boundary"]["fixed_point_materialized"] is False
    assert "mixed2" in layer["capability"]
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
