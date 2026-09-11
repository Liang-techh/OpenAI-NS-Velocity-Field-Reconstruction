from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_cutoff_schedule import (
    build_slow_borel_cutoff_schedule,
)
from openai_ns_reconstruction.background_cutoff_support import (
    classify_cutoff_support,
    exact_plateau_edge,
    exact_zero_edge,
)
from openai_ns_reconstruction.cutoffs import standard_cutoff


def _doubling_schedule():
    # For h=1 and C[j,m]=1, every local positive-order witness is 2 and
    # the recursive envelope is therefore exactly (1, 2, 4, 8).
    return build_slow_borel_cutoff_schedule(
        1.0,
        [
            [1.0] * 4,
            [1.0] * 5,
            [1.0] * 6,
        ],
    )


def test_classification_matches_independent_cutoff_values():
    schedule = _doubling_schedule()
    assert schedule.scales == (1, 2, 4, 8)

    cert = classify_cutoff_support(schedule, 0.2)
    assert cert.plateau_through == 1
    assert cert.transition_order == 2
    assert cert.zero_from == 3
    assert cert.stable_truncation_order_within_prefix == 2
    assert [cert.status(j) for j in range(4)] == [
        "leading",
        "plateau",
        "transition",
        "zero",
    ]

    # Independent executable cutoff path, not the certificate implementation.
    weights = [standard_cutoff(schedule.scales[j] * 0.2) for j in range(1, 4)]
    assert weights[0] == 1.0
    assert 0.0 < weights[1] < 1.0
    assert weights[2] == 0.0

    assert not cert.certifies_prefix_stability(1)
    assert cert.certifies_prefix_stability(2)
    assert cert.certifies_prefix_stability(3)


def test_closed_plateau_and_zero_boundaries_are_exact():
    schedule = _doubling_schedule()

    # q=1/4 gives a_1 q=1/2 and a_2 q=1 exactly.  Both numbers are
    # binary64-exact, so the Fraction-based classifier has no tolerance issue.
    cert = classify_cutoff_support(schedule, 0.25)
    assert cert.plateau_through == 1
    assert cert.transition_order is None
    assert cert.zero_from == 2
    assert cert.status(1) == "plateau"
    assert cert.status(2) == "zero"
    assert standard_cutoff(2.0 * 0.25) == 1.0
    assert standard_cutoff(4.0 * 0.25) == 0.0


def test_finite_prefix_ends_in_transition_fails_closed_on_tail():
    schedule = _doubling_schedule()
    cert = classify_cutoff_support(schedule, 0.1)

    assert cert.plateau_through == 2
    assert cert.transition_order == 3
    assert cert.zero_from is None
    assert not cert.tail_zero_certified_within_prefix
    assert cert.stable_truncation_order_within_prefix is None
    assert cert.status(3) == "transition"
    assert all(not cert.certifies_prefix_stability(j) for j in range(4))


def test_all_positive_orders_can_be_forced_zero_immediately():
    schedule = _doubling_schedule()
    cert = classify_cutoff_support(schedule, 0.6)

    assert cert.plateau_through == 0
    assert cert.transition_order is None
    assert cert.zero_from == 1
    assert cert.stable_truncation_order_within_prefix == 0
    assert cert.status(0) == "leading"
    assert all(cert.status(j) == "zero" for j in range(1, 4))


def test_exact_schedule_edges_are_rational():
    schedule = _doubling_schedule()
    assert exact_plateau_edge(schedule, 1) == Fraction(1, 4)
    assert exact_zero_edge(schedule, 1) == Fraction(1, 2)
    assert exact_plateau_edge(schedule, 3) == Fraction(1, 16)
    assert exact_zero_edge(schedule, 3) == Fraction(1, 8)


def test_fail_closed_inputs():
    schedule = _doubling_schedule()

    for q in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(ValueError):
            classify_cutoff_support(schedule, q)

    with pytest.raises(TypeError):
        classify_cutoff_support(object(), 0.25)
    with pytest.raises(ValueError):
        exact_plateau_edge(schedule, 0)
    with pytest.raises(ValueError):
        exact_zero_edge(schedule, schedule.max_order + 1)

    cert = classify_cutoff_support(schedule, 0.2)
    with pytest.raises(ValueError):
        cert.status(schedule.max_order + 1)
    with pytest.raises(ValueError):
        cert.certifies_prefix_stability(-1)
