import math

import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_limit_majorant import EndpointPowerLawMajorant


def _spatial_amplitude(x: float, y: float, z: float) -> float:
    return 1.0 + x * x + y * y + z * z


def _manufactured_jet(t: float, x: float, y: float, z: float, beta: float) -> np.ndarray:
    limit = np.array([x - y, z, x + y + z])
    correction = np.array([_spatial_amplitude(x, y, z), 0.0, 0.0])
    return limit + math.pow(1.0 - t, beta) * correction


def _manufactured_limit(x: float, y: float, z: float) -> np.ndarray:
    return np.array([x - y, z, x + y + z])


def test_power_majorant_controls_manufactured_uniform_endpoint_limit() -> None:
    beta = 0.4
    alpha = 1.0 - beta
    amplitude_sup = 1.75  # exact max of 1+x^2+y^2+z^2 on [-1/2,1/2]^3
    majorant = EndpointPowerLawMajorant(
        derivative_degree=3,
        spatial_window=2,
        coefficient=beta * amplitude_sup,
        singularity_exponent=alpha,
        valid_from=0.75,
        endpoint=1.0,
    )

    t0, t1 = 0.81, 0.973
    certificate = majorant.cauchy_certificate(t0, t1)
    assert certificate.certified

    actual_interval_sup = 0.0
    actual_tail_sup = 0.0
    for x in np.linspace(-0.5, 0.5, 5):
        for y in np.linspace(-0.5, 0.5, 5):
            for z in np.linspace(-0.5, 0.5, 5):
                before = _manufactured_jet(t0, x, y, z, beta)
                after = _manufactured_jet(t1, x, y, z, beta)
                limit = _manufactured_limit(x, y, z)
                actual_interval_sup = max(
                    actual_interval_sup, float(np.linalg.norm(after - before))
                )
                actual_tail_sup = max(
                    actual_tail_sup, float(np.linalg.norm(before - limit))
                )

    assert actual_interval_sup <= certificate.interval_upper_bound * (1.0 + 2e-14)
    assert actual_tail_sup <= certificate.endpoint_tail_upper_bound * (1.0 + 2e-14)

    # The cube corners attain the supplied amplitude bound, so the manufactured
    # example is a sharp cross-check of the antiderivative arithmetic rather
    # than merely a very loose inequality.
    expected_tail = amplitude_sup * math.pow(1.0 - t0, beta)
    assert math.isclose(
        certificate.endpoint_tail_upper_bound,
        expected_tail,
        rel_tol=2e-15,
        abs_tol=0.0,
    )


def test_interval_formula_matches_independent_midpoint_integration() -> None:
    majorant = EndpointPowerLawMajorant(
        derivative_degree=1,
        spatial_window=0,
        coefficient=2.3,
        singularity_exponent=0.35,
        valid_from=0.75,
    )
    t0, t1 = 0.8, 0.96
    analytic = majorant.interval_bound(t0, t1)

    # Independent numerical integration of C(1-t)^(-alpha); this is only a
    # regression oracle for the arithmetic implementation, not a proof that an
    # actual residual satisfies the supplied bound.
    count = 20000
    width = (t1 - t0) / count
    midpoints = t0 + (np.arange(count, dtype=float) + 0.5) * width
    numeric = float(
        np.sum(
            majorant.coefficient
            * np.power(majorant.endpoint - midpoints, -majorant.singularity_exponent)
        )
        * width
    )
    assert math.isclose(analytic, numeric, rel_tol=2e-9, abs_tol=1e-12)


def test_epsilon_entry_time_is_a_fail_closed_modulus_not_a_fitted_rate() -> None:
    majorant = EndpointPowerLawMajorant(
        derivative_degree=4,
        spatial_window=3,
        coefficient=0.8,
        singularity_exponent=0.5,
        valid_from=0.75,
    )
    epsilon = 0.12
    entry = majorant.epsilon_entry_time(epsilon)
    assert majorant.valid_from <= entry < majorant.endpoint
    assert majorant.endpoint_tail_bound(entry) <= epsilon * (1.0 + 4e-15)

    previous = entry - 1e-6
    if previous >= majorant.valid_from:
        assert majorant.endpoint_tail_bound(previous) > epsilon

    zero = EndpointPowerLawMajorant(
        derivative_degree=0,
        spatial_window=0,
        coefficient=0.0,
        singularity_exponent=0.9,
        valid_from=0.75,
    )
    assert zero.epsilon_entry_time(1e-12) == zero.valid_from
    assert zero.endpoint_tail_bound(0.9) == 0.0


def test_nonintegrable_or_out_of_scope_inputs_fail_closed() -> None:
    with pytest.raises(ValueError, match="0 <= alpha < 1"):
        EndpointPowerLawMajorant(0, 0, 1.0, 1.0)
    with pytest.raises(ValueError, match="0 <= alpha < 1"):
        EndpointPowerLawMajorant(0, 0, 1.0, -0.1)
    with pytest.raises(ValueError, match="coefficient must be nonnegative"):
        EndpointPowerLawMajorant(0, 0, -1.0, 0.5)
    with pytest.raises(ValueError, match="nonnegative integer"):
        EndpointPowerLawMajorant(True, 0, 1.0, 0.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="strictly before endpoint"):
        EndpointPowerLawMajorant(0, 0, 1.0, 0.5, valid_from=1.0, endpoint=1.0)

    majorant = EndpointPowerLawMajorant(0, 0, 1.0, 0.5, valid_from=0.75)
    with pytest.raises(ValueError, match="valid_from"):
        majorant.endpoint_tail_bound(0.7)
    with pytest.raises(ValueError, match="valid_from"):
        majorant.endpoint_tail_bound(1.0)
    with pytest.raises(ValueError, match="t0 <= t1"):
        majorant.interval_bound(0.9, 0.8)
    with pytest.raises(ValueError, match="epsilon must be positive"):
        majorant.epsilon_entry_time(0.0)
