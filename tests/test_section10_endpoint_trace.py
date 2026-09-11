import math

import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_limit_majorant import EndpointPowerLawMajorant
from openai_ns_reconstruction.section10_endpoint_ladder import Section10EndpointMajorantLadder
from openai_ns_reconstruction.section10_endpoint_trace import (
    section10_endpoint_trace_enclosures,
)


def _analytic_family():
    # Independent closed-form model: J_n(t) = L_n + a_n (1-t)^beta_n e_n.
    # Its time derivative has exact Euclidean norm
    # beta_n |a_n| (1-t)^(-(1-beta_n)), so the integrated majorant radius is
    # exactly |a_n| (1-t)^beta_n.
    betas = (1.0, 0.75, 0.5)
    amplitudes = (0.4, -0.7, 1.1)
    limits = []
    directions = []
    for degree in range(3):
        shape = (3,) + (4,) * degree
        limit = np.zeros(shape)
        limit.reshape(-1)[0] = 0.25 * (degree + 1)
        direction = np.zeros(shape)
        direction.reshape(-1)[-1] = 1.0
        limits.append(limit)
        directions.append(direction)

    def past(degree, x, y, z, t):
        del x, y, z
        return (
            limits[degree]
            + amplitudes[degree]
            * math.pow(1.0 - t, betas[degree])
            * directions[degree]
        )

    # Use a strict analytic upper bound (1% slack) rather than relying on
    # floating-point equality at the boundary of the closed ball.
    majorants = [
        EndpointPowerLawMajorant(
            derivative_degree=degree,
            spatial_window=3,
            coefficient=1.01 * betas[degree] * abs(amplitudes[degree]),
            singularity_exponent=1.0 - betas[degree],
            valid_from=0.75,
            endpoint=1.0,
        )
        for degree in range(3)
    ]
    return past, tuple(limits), tuple(amplitudes), tuple(betas), majorants


def test_endpoint_trace_enclosures_contain_independently_known_limits() -> None:
    past, limits, amplitudes, betas, majorants = _analytic_family()
    ladder = Section10EndpointMajorantLadder(majorants)
    t = 0.91
    enclosures = section10_endpoint_trace_enclosures(ladder, past, 0.1, -0.2, 0.3, t)

    assert tuple(item.degree for item in enclosures) == (0, 1, 2)
    assert all(item.certified for item in enclosures)
    for degree, enclosure in enumerate(enclosures):
        exact_remainder = abs(amplitudes[degree]) * math.pow(1.0 - t, betas[degree])
        expected_radius = 1.01 * exact_remainder
        assert np.isclose(enclosure.error_radius, expected_radius, rtol=3e-15, atol=0.0)
        # The endpoint tensor is known from the closed form, not produced by
        # the enclosure implementation or the majorant arithmetic.
        assert enclosure.contains(limits[degree])
        assert np.isclose(
            enclosure.distance_to(limits[degree]), exact_remainder, rtol=3e-15, atol=1e-16
        )


def test_endpoint_trace_enclosures_shrink_toward_t_one() -> None:
    past, _, _, _, majorants = _analytic_family()
    ladder = Section10EndpointMajorantLadder(majorants)
    early = section10_endpoint_trace_enclosures(ladder, past, 0.0, 0.0, 0.0, 0.80)
    late = section10_endpoint_trace_enclosures(ladder, past, 0.0, 0.0, 0.0, 0.98)
    assert all(b.error_radius < a.error_radius for a, b in zip(early, late))


def test_endpoint_trace_candidate_mismatch_fails_membership_without_relabeling_limit() -> None:
    past, limits, _, _, majorants = _analytic_family()
    ladder = Section10EndpointMajorantLadder(majorants)
    enclosure = section10_endpoint_trace_enclosures(
        ladder, past, 0.0, 0.0, 0.0, 0.95
    )[1]
    wrong = limits[1].copy()
    wrong.reshape(-1)[1] = 10.0
    assert not enclosure.contains(wrong)


def test_endpoint_trace_fails_closed_on_invalid_time_or_tensor_shape() -> None:
    past, _, _, _, majorants = _analytic_family()
    ladder = Section10EndpointMajorantLadder(majorants)
    with pytest.raises(ValueError, match="common_valid_from"):
        section10_endpoint_trace_enclosures(ladder, past, 0.0, 0.0, 0.0, 0.70)

    def bad_shape(degree, x, y, z, t):
        del degree, x, y, z, t
        return np.zeros(3)

    with pytest.raises(ValueError, match="order-1 full spacetime jet"):
        section10_endpoint_trace_enclosures(ladder, bad_shape, 0.0, 0.0, 0.0, 0.90)
