import math

import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_limit_majorant import EndpointPowerLawMajorant
from openai_ns_reconstruction.section10_endpoint_ladder import Section10EndpointMajorantLadder
from openai_ns_reconstruction.section10_endpoint_trace_consistency import (
    certify_section10_endpoint_trace_consistency,
)


def _ladder() -> Section10EndpointMajorantLadder:
    # J_n(t) below has amplitude c_n sqrt(1-t), hence
    # ||partial_t J_n|| = (c_n/2) (1-t)^(-1/2).
    amplitudes = (1.0, 0.5, 0.25)
    return Section10EndpointMajorantLadder(
        EndpointPowerLawMajorant(
            derivative_degree=n,
            spatial_window=2,
            coefficient=0.5 * amplitude,
            singularity_exponent=0.5,
            valid_from=0.75,
            endpoint=1.0,
        )
        for n, amplitude in enumerate(amplitudes)
    )


def _analytic_jets(degree: int, x: float, y: float, z: float, t: float) -> np.ndarray:
    amplitudes = (1.0, 0.5, 0.25)
    shape = (3,) + (4,) * degree
    out = np.zeros(shape)
    # A time-independent base verifies that only the prescribed endpoint tail
    # contributes to the difference between the two samples.
    out[(0,) + (0,) * degree] = 2.0 + x - 0.5 * y + z
    out[(1,) + (0,) * degree] = amplitudes[degree] * math.sqrt(1.0 - t)
    return out


def test_consistency_matches_independent_closed_form_ftc_budget() -> None:
    ladder = _ladder()
    t0, t1 = 0.80, 0.95
    certificate = certify_section10_endpoint_trace_consistency(
        ladder, _analytic_jets, 0.1, -0.2, 0.3, t0, t1
    )

    assert certificate.formal_consistency_ready
    assert [row.degree for row in certificate.rows] == [0, 1, 2]
    for degree, row in enumerate(certificate.rows):
        amplitude = (1.0, 0.5, 0.25)[degree]
        exact = amplitude * (math.sqrt(1.0 - t0) - math.sqrt(1.0 - t1))
        assert row.observed_distance == pytest.approx(exact, rel=2e-15, abs=2e-15)
        assert row.cauchy_upper_bound == pytest.approx(exact, rel=2e-15, abs=2e-15)
        assert row.compatible

    assert not certificate.actual_residual_jet_family_verified
    assert not certificate.uniform_spatial_majorant_verified
    assert not certificate.locally_uniform_endpoint_limits_verified
    assert not certificate.all_order_borel_smoothness_verified
    assert not certificate.paper_exact_velocity_available


def test_inconsistent_samples_fail_instead_of_refitting_the_majorant() -> None:
    ladder = _ladder()

    def bad_jets(degree: int, x: float, y: float, z: float, t: float) -> np.ndarray:
        out = _analytic_jets(degree, x, y, z, t).copy()
        if degree == 1 and t > 0.9:
            out[(2,) + (0,) * degree] = 1.0
        return out

    with pytest.raises(ValueError, match="derivative degree 1"):
        certify_section10_endpoint_trace_consistency(
            ladder, bad_jets, 0.0, 0.0, 0.0, 0.80, 0.95
        )


def test_wrong_jet_shape_and_time_window_fail_closed() -> None:
    ladder = _ladder()

    def wrong_shape(degree: int, x: float, y: float, z: float, t: float) -> np.ndarray:
        if degree == 2:
            return np.zeros((3, 4))
        return _analytic_jets(degree, x, y, z, t)

    with pytest.raises(ValueError, match="order-2 full spacetime jet"):
        certify_section10_endpoint_trace_consistency(
            ladder, wrong_shape, 0.0, 0.0, 0.0, 0.80, 0.90
        )

    with pytest.raises(ValueError, match="common_valid_from"):
        certify_section10_endpoint_trace_consistency(
            ladder, _analytic_jets, 0.0, 0.0, 0.0, 0.70, 0.90
        )
