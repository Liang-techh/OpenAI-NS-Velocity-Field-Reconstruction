from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_scale_schedule import (
    SpatialBorelScaleSchedule,
    borel_extension_from_template_bounds,
)


def _manual_bound_sum(provider, degree: int) -> Fraction:
    return sum(
        (Fraction(provider(m, degree, k)) for m in range(degree) for k in range(degree)),
        Fraction(0),
    )


def test_local_scale_is_exact_strict_witness_for_pinned_bound_sum() -> None:
    def bounds(m: int, j: int, k: int) -> Fraction:
        return Fraction((m + 1) * (k + 2), j + 3)

    schedule = SpatialBorelScaleSchedule(bounds)
    previous_scale = None
    for j in range(7):
        expected_sum = _manual_bound_sum(bounds, j)
        target = (1 << j) * expected_sum
        expected_local = target.numerator // target.denominator + 1

        cert = schedule.degree_certificate(j)
        assert cert.bound_sum == expected_sum
        assert cert.local_scale == expected_local
        assert cert.local_scale > target
        assert cert.scale >= cert.local_scale
        if previous_scale is not None:
            assert cert.scale >= 2 * previous_scale
        assert cert.bound_sum <= Fraction(cert.scale, 1 << j)
        assert cert.certified
        previous_scale = cert.scale


def test_tail_certificate_recomputes_localized_derivative_bound_exactly() -> None:
    def bounds(m: int, j: int, k: int) -> Fraction:
        return Fraction(m + k + 1, 5)

    schedule = SpatialBorelScaleSchedule(bounds)
    for j in range(1, 7):
        for m in range(j):
            for k in range(j):
                cert = schedule.tail_certificate(m, j, k)
                independent = Fraction(bounds(m, j, k), schedule.scale(j) ** (j - k))
                assert cert.direct_upper_bound == independent
                assert cert.geometric_target == Fraction(1, 1 << j)
                assert independent <= cert.geometric_target
                assert cert.certified


def test_schedule_is_wired_to_existing_full_jet_borel_path() -> None:
    def bounds(m: int, j: int, k: int) -> Fraction:
        return Fraction(1 + m + k, 8)

    def full_jets(j: int, x: float, y: float, z: float):
        tensor = np.zeros((3,) + (4,) * j)
        index = (slice(None),) + (0,) * j
        tensor[index] = np.array([j + 1.0, x - y, z])
        return tensor

    extension, schedule = borel_extension_from_template_bounds(full_jets, bounds)
    for j in range(6):
        # The extension applies its own doubling envelope to the same deterministic
        # localScale witness, so the two scale paths must agree exactly.
        assert extension.scale(j) == schedule.scale(j)

    value = extension(0.2, -0.1, 0.3, 1.0)
    assert np.array_equal(value, np.array([1.0, 0.3, 0.3]))


def test_bound_provider_is_cached_and_only_finite_triangle_is_needed_per_degree() -> None:
    calls: list[tuple[int, int, int]] = []

    def bounds(m: int, j: int, k: int) -> int:
        calls.append((m, j, k))
        return 1

    schedule = SpatialBorelScaleSchedule(bounds)
    assert schedule.bound_sum(3) == 9
    assert schedule.bound_sum(3) == 9
    assert sorted(calls) == sorted((m, 3, k) for m in range(3) for k in range(3))


def test_invalid_template_bounds_and_tail_indices_fail_closed() -> None:
    with pytest.raises(ValueError, match="callable"):
        SpatialBorelScaleSchedule(None)  # type: ignore[arg-type]

    negative = SpatialBorelScaleSchedule(lambda m, j, k: -1)
    with pytest.raises(ValueError, match="finite nonnegative"):
        negative.bound_sum(2)

    nonfinite = SpatialBorelScaleSchedule(lambda m, j, k: np.nan)
    with pytest.raises(ValueError, match="finite nonnegative"):
        nonfinite.bound_sum(1)

    schedule = SpatialBorelScaleSchedule(lambda m, j, k: 0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        schedule.bound_sum(True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="requires spatial window < degree"):
        schedule.tail_certificate(1, 1, 0)
    with pytest.raises(ValueError, match="requires spatial window < degree"):
        schedule.tail_certificate(0, 1, 1)
