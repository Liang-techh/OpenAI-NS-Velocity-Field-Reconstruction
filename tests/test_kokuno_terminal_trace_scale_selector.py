from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_terminal_trace_scale_selector import (
    TerminalTraceScaleError,
    cutoff_derivative_constant,
    initial_terminal_scale,
    select_terminal_scale,
    terminal_constraints,
    terminal_monomial_bound,
)


CHI = {0: Fraction(1), 1: Fraction(2), 2: Fraction(3)}
JETS = {0: Fraction(1), 1: Fraction(5, 2), 2: Fraction(7, 3)}


def test_initial_scale_is_exact_public_b0():
    assert initial_terminal_scale() == 1


def test_cutoff_constants_match_independent_finite_sum():
    assert cutoff_derivative_constant(
        stage_index=4, time_order=0, cutoff_derivative_norms=CHI
    ) == Fraction(1, 24)
    assert cutoff_derivative_constant(
        stage_index=4, time_order=1, cutoff_derivative_norms=CHI
    ) == Fraction(1, 4)
    assert cutoff_derivative_constant(
        stage_index=4, time_order=2, cutoff_derivative_norms=CHI
    ) == Fraction(31, 24)


def test_complete_constraint_family_has_expected_triangular_size():
    constraints = terminal_constraints(
        stage_index=4,
        scale=5,
        cutoff_derivative_norms=CHI,
        jet_norms=JETS,
    )
    assert len(constraints) == 6
    assert {(item.spatial_order, item.time_order) for item in constraints} == {
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (2, 0),
    }


def test_selector_returns_least_increasing_integer_and_exact_margin():
    selection = select_terminal_scale(
        stage_index=4,
        previous_scale=2,
        cutoff_derivative_norms=CHI,
        jet_norms=JETS,
        max_scale=20,
    )
    assert selection.scale == 5
    assert all(item.satisfied for item in selection.constraints)
    assert selection.max_target_ratio == Fraction(62, 75)

    predecessor = terminal_constraints(
        stage_index=4,
        scale=4,
        cutoff_derivative_norms=CHI,
        jet_norms=JETS,
    )
    assert any(not item.satisfied for item in predecessor)
    worst = max(predecessor, key=lambda item: item.bound / item.target)
    assert (worst.spatial_order, worst.time_order) == (0, 2)
    assert worst.bound == Fraction(31, 384)
    assert worst.target == Fraction(1, 16)


def test_single_constraint_exposes_exact_source_power_and_target():
    item = terminal_monomial_bound(
        stage_index=4,
        spatial_order=1,
        time_order=1,
        scale=5,
        cutoff_derivative_norms=CHI,
        jet_norms=JETS,
    )
    assert item.cutoff_constant == Fraction(1, 4)
    assert item.jet_norm == Fraction(5, 2)
    assert item.bound == Fraction(1, 200)
    assert item.target == Fraction(1, 16)
    assert item.satisfied


def test_exact_small_mutation_is_not_hidden_by_tolerance():
    item = terminal_monomial_bound(
        stage_index=4,
        spatial_order=0,
        time_order=2,
        scale=5,
        cutoff_derivative_norms=CHI,
        jet_norms=JETS,
    )
    mutated_jets = dict(JETS)
    mutated_jets[0] += Fraction(1, 2**40)
    mutated = terminal_monomial_bound(
        stage_index=4,
        spatial_order=0,
        time_order=2,
        scale=5,
        cutoff_derivative_norms=CHI,
        jet_norms=mutated_jets,
    )
    assert mutated.bound - item.bound == Fraction(31, 24 * 25 * 2**40)
    assert mutated.bound != item.bound


def test_fail_closed_on_inexact_inputs_missing_bounds_and_search_exhaustion():
    with pytest.raises(TypeError):
        select_terminal_scale(
            stage_index=4,
            previous_scale=2,
            cutoff_derivative_norms={0: Fraction(1), 1: Fraction(2), 2: 3.0},
            jet_norms=JETS,
            max_scale=20,
        )
    with pytest.raises(TerminalTraceScaleError):
        select_terminal_scale(
            stage_index=4,
            previous_scale=2,
            cutoff_derivative_norms={0: Fraction(1), 1: Fraction(2)},
            jet_norms=JETS,
            max_scale=20,
        )
    with pytest.raises(TerminalTraceScaleError):
        select_terminal_scale(
            stage_index=4,
            previous_scale=2,
            cutoff_derivative_norms=CHI,
            jet_norms=JETS,
            max_scale=4,
        )
    with pytest.raises(TypeError):
        select_terminal_scale(
            stage_index=True,
            previous_scale=2,
            cutoff_derivative_norms=CHI,
            jet_norms=JETS,
            max_scale=20,
        )
