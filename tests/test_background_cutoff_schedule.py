import math
import sys

import pytest

from openai_ns_reconstruction.background_cutoff_schedule import (
    SlowBorelCutoffSchedule,
    build_slow_borel_cutoff_schedule,
    build_slow_borel_cutoff_schedule_from_provider,
    local_scale_from_template_bounds,
)


def test_local_scales_match_independent_closed_form_cases() -> None:
    # For p=0 the Lean smallness row is C * b^(-h*j) <= 2^(-j).
    # These choices make the exact real thresholds 16, 4 and 16.
    assert local_scale_from_template_bounds(0.5, 1, [2.0] * 4) == 16
    assert local_scale_from_template_bounds(0.5, 2, [1.0] * 5) >= 4
    assert local_scale_from_template_bounds(0.5, 3, [8.0] * 6) >= 16


def test_doubling_envelope_and_all_edge_majorants() -> None:
    rows = [
        [2.0, 1.0, 0.5, 2.0],
        [1.0, 0.75, 0.5, 0.25, 1.0],
        [8.0, 1.0, 2.0, 4.0, 8.0, 0.5],
    ]
    schedule = build_slow_borel_cutoff_schedule(
        0.5, rows, initial_lower_bound=3
    )

    # The local positive-order witnesses are 16, at least 4, and at least 16;
    # the paper/Lean doubling envelope then forces 3,16,32,64.
    assert schedule.scales == (3, 16, 32, 64)
    assert schedule.local_scales[1] == 16
    assert schedule.local_scales[2] <= 32
    assert schedule.local_scales[3] <= 64

    for j, row in enumerate(rows, start=1):
        assert schedule.scales[j] >= 2 * schedule.scales[j - 1]
        assert schedule.scales[j] >= schedule.local_scales[j]
        edge = 1.0 / schedule.scales[j]
        for m, c in enumerate(row):
            # Independent direct evaluation of the SlowBorelBase p=0 edge row.
            lhs = c * edge ** (0.5 * j)
            rhs = 0.5 ** j
            assert lhs <= rhs * (1.0 + 2e-15)
            assert schedule.edge_log_margin(j, m) >= -1e-14
            got_lhs, got_rhs = schedule.majorant_at(j, m, 0.37 * edge)
            assert math.isclose(got_lhs, c * (0.37 * edge) ** (0.5 * j), rel_tol=2e-15)
            assert got_rhs == rhs
            assert got_lhs <= got_rhs


def test_wide_integer_scale_crosses_binary64_without_clipping() -> None:
    # This is a representation-layer regression, not paper coefficient data.
    # C=1e308, h=1/2 and j=1 require log(b) >= 2(log(2)+log(C)),
    # far beyond log(max_float), while Python's exact integer scale remains valid.
    row = [1.0e308] * 4
    scale = local_scale_from_template_bounds(0.5, 1, row)
    assert scale > sys.float_info.max

    lhs_log = math.log(row[0]) - 0.5 * math.log(scale)
    assert lhs_log <= -math.log(2.0) + 32 * sys.float_info.epsilon

    schedule = build_slow_borel_cutoff_schedule(0.5, [row])
    assert schedule.scales[1] == scale
    assert schedule.edge_log_margin(1, 0) >= -32 * sys.float_info.epsilon
    assert schedule.reciprocal_support_log_edge(1) == -math.log(scale)

    # 1/scale is below the least positive binary64 subnormal.  The executable
    # accessor must not silently return 0 and pretend that is an exact edge.
    with pytest.raises(OverflowError, match="underflows binary64"):
        schedule.reciprocal_support_edge(1)


def test_provider_is_called_for_exact_triangular_jet_set() -> None:
    calls: list[tuple[int, int]] = []

    def provider(j: int, m: int) -> float:
        calls.append((j, m))
        return 1.0 + j + m / 10.0

    schedule = build_slow_borel_cutoff_schedule_from_provider(
        0.2, 3, provider, initial_lower_bound=5
    )
    assert schedule.max_order == 3
    assert calls == [
        (1, 0), (1, 1), (1, 2), (1, 3),
        (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
        (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
    ]
    assert schedule.scales[0] >= 5
    assert all(schedule.scales[j] >= 2 * schedule.scales[j - 1] for j in range(1, 4))


def test_schedule_rejects_nonanalytic_or_malformed_inputs() -> None:
    with pytest.raises(ValueError, match="order zero"):
        local_scale_from_template_bounds(0.1, 0, [1.0, 1.0, 1.0])
    with pytest.raises(ValueError, match=r"exactly j\+3"):
        local_scale_from_template_bounds(0.1, 1, [1.0, 1.0, 1.0])
    with pytest.raises(ValueError, match="finite and positive"):
        local_scale_from_template_bounds(0.1, 1, [1.0, 1.0, 1.0, 0.0])
    with pytest.raises(ValueError, match="finite and positive"):
        build_slow_borel_cutoff_schedule(0.0, [])
    with pytest.raises(ValueError, match="nonnegative integer"):
        build_slow_borel_cutoff_schedule(0.1, [], initial_lower_bound=-1)
    with pytest.raises(TypeError, match="callable"):
        build_slow_borel_cutoff_schedule_from_provider(0.1, 2, 3.0)  # type: ignore[arg-type]


def test_dataclass_fails_closed_if_envelope_or_majorant_is_forged() -> None:
    rows = ((2.0, 2.0, 2.0, 2.0),)
    with pytest.raises(ValueError, match="at least double"):
        SlowBorelCutoffSchedule(0.5, 3, rows, (3, 1), (3, 5))
    with pytest.raises(ValueError, match="dyadic edge bound"):
        SlowBorelCutoffSchedule(0.5, 1, rows, (1, 1), (1, 2))


def test_majorant_domain_is_fail_closed() -> None:
    schedule = build_slow_borel_cutoff_schedule(0.5, [[2.0] * 4])
    edge = schedule.reciprocal_support_edge(1)
    with pytest.raises(ValueError, match="outside"):
        schedule.majorant_at(1, 0, math.nextafter(edge, math.inf))
    with pytest.raises(ValueError, match=r"m <= j\+2"):
        schedule.majorant_at(1, 4, 0.5 * edge)
    with pytest.raises(ValueError, match="positive orders"):
        schedule.edge_log_margin(0, 0)
