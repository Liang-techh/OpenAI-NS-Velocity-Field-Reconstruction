import math

import pytest
from scipy.integrate import solve_ivp

from openai_ns_reconstruction.outgoing_tail import (
    OutgoingCoreParameters,
    TailData,
    clock_weight,
    final_angular,
    power_constant,
    radial_amplitude,
    release_slope,
    shape_exponent,
    tail_shape,
)


def _data() -> TailData:
    return TailData(OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0), h=0.01)


def test_parameter_geometry_and_fail_closed() -> None:
    d = _data()
    assert d.core.drop_length == pytest.approx(math.e + 10.0)
    assert d.core.pulse_length == pytest.approx(260.0)
    assert d.flatten_end > d.core.endpoint
    assert d.release_start > d.flatten_end
    assert d.ramp_end > d.second_ramp_start
    assert 0.0 < d.rho < 0.5
    with pytest.raises(ValueError):
        OutgoingCoreParameters(P=0.0, m=1.0, lam=0.05, wait=30.0)
    with pytest.raises(ValueError):
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.1, wait=30.0)
    with pytest.raises(ValueError):
        TailData(d.core, h=0.03)


def test_final_angular_matches_exact_ideal_prefix() -> None:
    d = _data()
    for y in (-8.0, -1.0, 0.0):
        for eta in (-2.0, 0.0, 1.5):
            expected = d.core.P * math.exp(y / 10.0) / (1.0 + eta * eta)
            assert final_angular(d, y, eta) == pytest.approx(expected, rel=3e-13, abs=1e-14)
        assert clock_weight(d, y) == pytest.approx(d.core.P**2 * math.exp(y / 5.0), rel=5e-13)
        assert shape_exponent(d, y) == 1.0


def test_flattened_uniform_wait_is_eta_independent() -> None:
    d = _data()
    y = d.flatten_end + 0.4 * d.uniform_wait
    assert y < d.release_start
    expected = radial_amplitude(d.core, y) / 2.0
    assert final_angular(d, y, 0.0) == pytest.approx(expected, rel=3e-12)
    assert final_angular(d, y, 3.0) == pytest.approx(expected, rel=3e-12)
    assert shape_exponent(d, y) == 0.0


def test_release_lag_endpoint_agrees_with_independent_ode_solver() -> None:
    d = _data()
    q0 = (d.core.lam - d.h) / (1.0 - d.core.lam)

    def rhs(t, q):
        slope = release_slope(d, float(t))
        return [-slope - d.h - (1.0 + slope) * q[0]]

    sol = solve_ivp(rhs, (0.0, d.ramp_end), [q0], rtol=2e-11, atol=2e-13, max_step=0.05)
    assert sol.success
    assert d.release_lag_at_ramp_end == pytest.approx(sol.y[0, -1], rel=2e-8, abs=2e-10)
    assert d.release_lag_at_ramp_end > d.tail_debt > 0.0
    assert d.release_lag_at_ramp_end * math.exp(-(1.0 - d.h) * d.decay_hold) == pytest.approx(
        d.tail_debt, rel=2e-13
    )


def test_terminal_taper_and_eventual_power_identity() -> None:
    d = _data()
    assert tail_shape(d, 1.0) == pytest.approx(1.0 - d.rho)
    assert tail_shape(d, 3.0) == pytest.approx(1.0)

    y1 = d.tail_end + 0.5
    y2 = y1 + 2.0
    expected1 = power_constant(d) * math.exp(-(0.5 + d.h) * y1)
    expected_ratio = math.exp(-(0.5 + d.h) * (y2 - y1))
    f1 = final_angular(d, y1, 0.0)
    f2 = final_angular(d, y2, 4.0)
    assert f1 == pytest.approx(expected1, rel=5e-11)
    assert f2 / f1 == pytest.approx(expected_ratio, rel=5e-11)
    assert final_angular(d, y1, 4.0) == pytest.approx(f1, rel=5e-12)
