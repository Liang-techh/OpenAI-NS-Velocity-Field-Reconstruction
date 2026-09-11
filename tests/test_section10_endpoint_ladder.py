import math

import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_limit_majorant import EndpointPowerLawMajorant
from openai_ns_reconstruction.section10_endpoint_ladder import Section10EndpointMajorantLadder
from openai_ns_reconstruction.time_localization import LATE_START, SECTION10_ENDPOINT


def _majorant(
    degree: int,
    coefficient: float,
    exponent: float,
    *,
    valid_from: float = LATE_START,
    window: int = 2,
    endpoint: float = SECTION10_ENDPOINT,
) -> EndpointPowerLawMajorant:
    return EndpointPowerLawMajorant(
        derivative_degree=degree,
        spatial_window=window,
        coefficient=coefficient,
        singularity_exponent=exponent,
        valid_from=valid_from,
        endpoint=endpoint,
    )


def test_ladder_builds_one_common_finite_order_endpoint_modulus() -> None:
    # Intentionally supply the entries out of order.  The ladder may reorder
    # them, but it must preserve exact contiguous degrees and use the latest
    # validity start as the common interval.
    ladder = Section10EndpointMajorantLadder(
        [
            _majorant(2, 1.5, 0.25, valid_from=0.80),
            _majorant(0, 0.5, 0.0, valid_from=0.75),
            _majorant(1, 2.0, 0.5, valid_from=0.78),
        ]
    )
    assert ladder.certified
    assert ladder.max_degree == 2
    assert ladder.spatial_window == 2
    assert ladder.endpoint == SECTION10_ENDPOINT
    assert ladder.common_valid_from == 0.80
    assert tuple(item.derivative_degree for item in ladder.majorants) == (0, 1, 2)

    t = 0.9
    # Independent closed-form oracle for C/(1-alpha)*(1-t)^(1-alpha).
    expected = tuple(
        item.coefficient
        / (1.0 - item.singularity_exponent)
        * math.pow(SECTION10_ENDPOINT - t, 1.0 - item.singularity_exponent)
        for item in ladder.majorants
    )
    actual = ladder.endpoint_tail_bounds(t)
    assert np.allclose(actual, expected, rtol=2e-15, atol=0.0)
    assert ladder.uniform_endpoint_tail_bound(t) == max(actual)

    later = 0.99
    assert ladder.uniform_endpoint_tail_bound(later) < ladder.uniform_endpoint_tail_bound(t)


def test_ladder_cauchy_bound_is_maximum_of_independent_antiderivatives() -> None:
    ladder = Section10EndpointMajorantLadder(
        [_majorant(0, 1.0, 0.0), _majorant(1, 3.0, 0.25)]
    )
    t0, t1 = 0.82, 0.94
    certificates = ladder.cauchy_certificates(t0, t1)
    assert all(item.certified for item in certificates)

    expected = []
    for item in ladder.majorants:
        power = 1.0 - item.singularity_exponent
        expected.append(
            item.coefficient
            / power
            * (
                math.pow(SECTION10_ENDPOINT - t0, power)
                - math.pow(SECTION10_ENDPOINT - t1, power)
            )
        )
    assert np.allclose(
        [item.interval_upper_bound for item in certificates],
        expected,
        rtol=2e-15,
        atol=0.0,
    )
    assert np.isclose(ladder.uniform_cauchy_bound(t0, t1), max(expected), rtol=2e-15)


def test_ladder_fails_closed_on_missing_or_duplicate_degrees() -> None:
    with pytest.raises(ValueError, match="0..N"):
        Section10EndpointMajorantLadder([_majorant(0, 1.0, 0.1), _majorant(2, 1.0, 0.1)])
    with pytest.raises(ValueError, match="0..N"):
        Section10EndpointMajorantLadder([_majorant(0, 1.0, 0.1), _majorant(0, 2.0, 0.2)])


def test_ladder_fails_closed_outside_official_endpoint_regime() -> None:
    with pytest.raises(ValueError, match="t>=3/4"):
        Section10EndpointMajorantLadder(
            [_majorant(0, 1.0, 0.2, valid_from=0.60)]
        )
    with pytest.raises(ValueError, match="endpoint=1"):
        Section10EndpointMajorantLadder(
            [_majorant(0, 1.0, 0.2, endpoint=1.2)]
        )
    with pytest.raises(ValueError, match="one spatial window"):
        Section10EndpointMajorantLadder(
            [_majorant(0, 1.0, 0.2, window=1), _majorant(1, 1.0, 0.2, window=2)]
        )


def test_ladder_rejects_times_before_common_validity() -> None:
    ladder = Section10EndpointMajorantLadder(
        [_majorant(0, 1.0, 0.2, valid_from=0.75), _majorant(1, 1.0, 0.2, valid_from=0.85)]
    )
    with pytest.raises(ValueError, match="common_valid_from"):
        ladder.uniform_endpoint_tail_bound(0.84)
    with pytest.raises(ValueError, match="t0 <= t1"):
        ladder.uniform_cauchy_bound(0.95, 0.90)
