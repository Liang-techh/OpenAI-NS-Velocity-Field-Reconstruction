import math

import numpy as np
import pytest

from openai_ns_reconstruction.primary_amplitude_ode import PrimaryModalDatum
from openai_ns_reconstruction.primary_amplitude_solution import (
    PrimaryAmplitudeInterval,
    primary_modal_rhs_at,
    solve_primary_amplitude_ivp,
    volterra_diagnostic,
)


def _constant_datum():
    # With a=b=c=0, h=1 and eigenRate=0, the reference-mode subtraction
    # makes the modal operator exactly -d I, d=j^2*viscosity.
    return PrimaryModalDatum(
        rho=0.0,
        rho_dot=0.0,
        g_K=0.0,
        F=0.0,
        N_theta=1.0,
        g_N=0.0,
        rotation=0.0,
        eigenvalue=1.7,
        eigenvector=1.0,
        eigen_rate=0.0,
        viscosity=0.2,
    )


def _forcing():
    # rho=0 gives forceX=-f_radial and forceY=-f_N.
    return (0.6, 0.0, -0.4)


def _exact_path(t, *, start, initial):
    d = 4.0 * 0.2
    # Modal forcing = ((-0.6+0.4)/2, (-0.6-0.4)/2) = (-0.1,-0.5).
    g = np.array([-0.1, -0.5])
    tau = t - start
    decay = math.exp(-d * tau)
    return decay * initial + (1.0 - decay) * g / d


def test_finite_interval_solver_matches_independent_constant_coefficient_oracle():
    interval = PrimaryAmplitudeInterval(0.25, 1.1, 2)
    initial = np.array([0.7, -0.3])
    data = _constant_datum()
    path = solve_primary_amplitude_ivp(
        interval=interval,
        initial_state=initial,
        datum_at=lambda _t: data,
        forcing_at=lambda _t: _forcing(),
        rtol=1e-12,
        atol=1e-14,
    )

    for t in (0.25, 0.4, 0.83, 1.1):
        expected = _exact_path(t, start=interval.start, initial=initial)
        assert np.allclose(path(t), expected, rtol=3e-11, atol=3e-12)


def test_volterra_diagnostic_accepts_exact_path_and_detects_perturbation():
    interval = PrimaryAmplitudeInterval(-0.2, 0.9, 2)
    initial = np.array([0.2, 0.9])
    data = _constant_datum()

    exact = lambda t: _exact_path(t, start=interval.start, initial=initial)
    diagnostic = volterra_diagnostic(
        interval=interval,
        initial_state=initial,
        state_at=exact,
        datum_at=lambda _t: data,
        forcing_at=lambda _t: _forcing(),
        time=0.73,
        epsabs=1e-12,
        epsrel=1e-12,
    )
    assert diagnostic.defect_norm < 2e-11

    def perturbed(t):
        value = exact(t).copy()
        value[0] += 1e-3 * (t - interval.start) ** 2
        return value

    bad = volterra_diagnostic(
        interval=interval,
        initial_state=initial,
        state_at=perturbed,
        datum_at=lambda _t: data,
        forcing_at=lambda _t: _forcing(),
        time=0.73,
        epsabs=1e-12,
        epsrel=1e-12,
    )
    assert bad.defect_norm > 1e-4


def test_rhs_and_interval_inputs_fail_closed():
    interval = PrimaryAmplitudeInterval(0.0, 1.0, 2)
    data = _constant_datum()

    with pytest.raises(ValueError, match="end"):
        PrimaryAmplitudeInterval(1.0, 0.0, 2)
    with pytest.raises(ValueError, match="harmonic"):
        PrimaryAmplitudeInterval(0.0, 1.0, True)

    with pytest.raises(TypeError, match="PrimaryModalDatum"):
        primary_modal_rhs_at(
            time=0.5,
            state=(0.0, 0.0),
            interval=interval,
            datum_at=lambda _t: object(),
            forcing_at=lambda _t: _forcing(),
        )
    with pytest.raises(ValueError, match="three finite scalars"):
        primary_modal_rhs_at(
            time=0.5,
            state=(0.0, 0.0),
            interval=interval,
            datum_at=lambda _t: data,
            forcing_at=lambda _t: (0.0, float("nan"), 0.0),
        )
    with pytest.raises(ValueError, match="closed amplitude interval"):
        primary_modal_rhs_at(
            time=1.1,
            state=(0.0, 0.0),
            interval=interval,
            datum_at=lambda _t: data,
            forcing_at=lambda _t: _forcing(),
        )


def test_volterra_checker_requires_exact_initial_value_not_fitted_endpoint():
    interval = PrimaryAmplitudeInterval(0.0, 0.5, 2)
    initial = np.array([0.2, 0.1])
    data = _constant_datum()

    with pytest.raises(ValueError, match=r"state_at\(start\)"):
        volterra_diagnostic(
            interval=interval,
            initial_state=initial,
            state_at=lambda t: initial + np.array([1e-12, 0.0]),
            datum_at=lambda _t: data,
            forcing_at=lambda _t: _forcing(),
            time=0.4,
        )
