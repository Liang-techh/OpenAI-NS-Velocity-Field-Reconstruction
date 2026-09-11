import math

import pytest

from openai_ns_reconstruction.natural_axis import A, D, real_phase
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_analytic_neighborhood import (
    certify_actual_schedule_analytic_inputs,
    certify_axis_analytic_neighborhood,
)
from openai_ns_reconstruction.schedule_axis_pressure import (
    axis_pressure,
    axis_pressure_derivative,
)


def _data(*, P: float = 2.0, h: float = 0.01) -> TailData:
    return TailData(
        OutgoingCoreParameters(P=P, m=1.0, lam=0.05, wait=30.0),
        h=h,
    )


def _complex_h(h: float, j: float, z: complex) -> complex:
    return D(h) * z + (1.0 - z * z) * (4.0 * z + j)


def test_actual_schedule_supplies_nonzero_tube_common_bound_and_phase_bound() -> None:
    cert = certify_actual_schedule_analytic_inputs(_data(), 0.05)
    n = cert.neighborhood

    assert n.radius > 0.0
    assert n.radius <= 0.25
    assert n.epsilon == n.radius / 2.0
    assert n.h_perturbation_fraction <= 0.25 * (1.0 + 1e-12)
    assert n.denominator_lower > 0.0
    assert n.common_field_sup_upper == pytest.approx(
        1.0 + math.fsum(n.field_bounds.values())
    )
    assert n.common_field_sup_upper > max(n.field_bounds.values())
    assert n.axis_phase_real_part_sup_upper > 0.0
    assert n.sigma == cert.low_z.cutoff.sigma
    assert n.clock_mass_upper == cert.low_z.clock_mass_upper


def test_complex_rational_fields_stay_inside_certified_bounds_on_independent_points() -> None:
    cert = certify_actual_schedule_analytic_inputs(_data(), 0.05)
    n = cert.neighborhood
    h = n.parameters.h
    j = n.parameters.j
    sigma = n.sigma

    # Sampling is only an independent regression consequence of the analytic
    # derivative/factor-separation certificate; it does not construct rho.
    directions = (1j, -1j, 0.6 + 0.8j, -0.6 + 0.8j)
    for x in (-1.1, -0.5, 0.0, 0.5, 1.1):
        for direction in directions:
            z = complex(x) + n.radius * direction
            hh = _complex_h(h, j, z)
            ll = 1.0 - 2.0 * h * z * z
            denominator = hh * hh + sigma * sigma
            chi = hh * hh / denominator
            gradient = -ll * hh / denominator

            assert abs(z.imag) < 0.5
            assert abs(ll) >= n.l_tube_lower * (1.0 - 1e-12)
            assert abs(denominator) >= n.denominator_lower * (1.0 - 1e-10)
            assert abs(chi) <= n.field_bounds["chi"] * (1.0 + 1e-10)
            assert abs(gradient) <= n.field_bounds["gradient"] * (1.0 + 1e-10)


def test_pressure_and_zstar_bounds_crosscheck_actual_real_schedule() -> None:
    data = _data()
    cert = certify_actual_schedule_analytic_inputs(data, 0.05)
    n = cert.neighborhood
    h = n.parameters.h
    j = n.parameters.j

    for eta in (-1.1, -0.5, 0.0, 0.5, 1.1):
        p = axis_pressure(data, eta)
        dp = axis_pressure_derivative(data, eta)
        z = float(eta)
        u = 4.0 * z + j
        hh = D(h) * z + (1.0 - z * z) * u
        zstar = (
            -A(h) * (1.0 - 2.0 * z * u) * u
            - 4.0 * hh
            - (1.0 - z * z) * dp
            + 4.0 * A(h) * z * p
        )

        assert abs(p) <= n.pressure_upper * (1.0 + 1e-10)
        assert abs(dp) <= n.pressure_derivative_upper * (1.0 + 1e-10)
        assert abs(zstar) <= n.field_bounds["zStar"] * (1.0 + 1e-10)

    for eta in (-1.0, -0.25, 0.25, 1.0):
        phase = real_phase(h, j, n.sigma, eta, samples=401)
        assert abs(phase) <= n.axis_phase_real_part_sup_upper * (1.0 + 1e-10)


def test_pure_certificate_fails_closed_outside_required_inputs() -> None:
    with pytest.raises(ValueError, match="1/100"):
        certify_axis_analytic_neighborhood(
            h=0.02,
            j=0.05,
            sigma=0.01,
            clock_mass_upper=1.0,
        )
    with pytest.raises(ValueError, match="1/20"):
        certify_axis_analytic_neighborhood(
            h=0.01,
            j=0.051,
            sigma=0.01,
            clock_mass_upper=1.0,
        )
    with pytest.raises(ValueError, match="sigma"):
        certify_axis_analytic_neighborhood(
            h=0.01,
            j=0.05,
            sigma=0.0,
            clock_mass_upper=1.0,
        )
    with pytest.raises(ValueError, match="clock_mass_upper"):
        certify_axis_analytic_neighborhood(
            h=0.01,
            j=0.05,
            sigma=0.01,
            clock_mass_upper=0.0,
        )
    with pytest.raises(ValueError, match="P >= 2"):
        certify_actual_schedule_analytic_inputs(_data(P=1.5), 0.05)
