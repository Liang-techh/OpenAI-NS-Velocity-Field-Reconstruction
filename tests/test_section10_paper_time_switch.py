from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section10_paper_time_switch import (
    ONE_PLATEAU_START,
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
    PINNED_SOURCE_FILE,
    Section10PaperTimeSwitchSource,
    ZERO_PLATEAU_ABS_RADIUS,
)


def test_pinned_time_switch_exact_plateau_geometry() -> None:
    source = Section10PaperTimeSwitchSource.pinned()

    assert source.source_key == (
        PINNED_LEAN_REPOSITORY,
        PINNED_LEAN_COMMIT,
        PINNED_SOURCE_FILE,
    )
    assert ZERO_PLATEAU_ABS_RADIUS == Fraction(3, 8)
    assert ONE_PLATEAU_START == Fraction(3, 4)

    assert source.certified_value(Fraction(-3, 8)) == 0
    assert source.certified_value(Fraction(0)) == 0
    assert source.certified_value(Fraction(3, 8)) == 0
    assert source.certified_value(Fraction(1, 2)) is None
    assert source.certified_value(Fraction(3, 4)) == 1
    assert source.certified_value(Fraction(1)) == 1


def test_time_switch_is_theorem_certified_inert_near_t1_without_claiming_extension() -> None:
    source = Section10PaperTimeSwitchSource.pinned()

    assert source.time_switch_contdiff_theorem_bound is True
    assert source.endpoint_t1_switch_value_certified_one is True
    assert source.endpoint_t1_switch_positive_derivatives_certified_zero is True
    assert source.positive_derivatives_certified_zero_late(Fraction(3, 4)) is False
    assert source.positive_derivatives_certified_zero_late(Fraction(3, 4) + Fraction(1, 1000)) is True

    assert source.section9_field_smooth_extension_through_t1_constructed is False
    assert source.endpoint_residual_closure_verified is False
    assert source.paper_exact_velocity_available is False
    assert source.actual_mathlib_bump_numerically_evaluated is False
    assert source.lean_theorems_machine_replayed_in_python is False


def test_derivative_support_collar_preserves_exact_nonnegative_boundary() -> None:
    source = Section10PaperTimeSwitchSource.pinned()

    assert source.positive_derivative_support_may_be_nonzero_on_nonnegative_axis(Fraction(0)) is False
    assert source.positive_derivative_support_may_be_nonzero_on_nonnegative_axis(Fraction(3, 8)) is True
    assert source.positive_derivative_support_may_be_nonzero_on_nonnegative_axis(Fraction(1, 2)) is True
    assert source.positive_derivative_support_may_be_nonzero_on_nonnegative_axis(Fraction(3, 4)) is True
    assert source.positive_derivative_support_may_be_nonzero_on_nonnegative_axis(Fraction(1)) is False

    with pytest.raises(ValueError, match="t >= 0"):
        source.positive_derivative_support_may_be_nonzero_on_nonnegative_axis(Fraction(-1, 10))


def test_arbitrary_time_window_or_source_drift_fails_closed() -> None:
    source = Section10PaperTimeSwitchSource.pinned()

    mutations = (
        {"lean_repository": "example/surrogate"},
        {"lean_commit": "0" * 40},
        {"source_file": "NavierStokes/OtherCutoffs.lean"},
        {"definition_name": "genericTimeWindow"},
        {"zero_theorem": "surrogate.zero"},
        {"one_theorem": "surrogate.one"},
        {"eventually_one_theorem": "surrogate.eventually_one"},
        {"late_derivative_theorem": "surrogate.derivative"},
        {"derivative_support_theorem": "surrogate.support"},
        {"zero_plateau_abs_radius": Fraction(1, 4)},
        {"one_plateau_start": Fraction(7, 8)},
    )
    for mutation in mutations:
        with pytest.raises((TypeError, ValueError)):
            replace(source, **mutation)


def test_exact_boundary_queries_reject_binary_float_shortcuts() -> None:
    source = Section10PaperTimeSwitchSource.pinned()

    with pytest.raises(TypeError, match="exact rational"):
        source.certified_value(0.75)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="exact rational"):
        source.positive_derivatives_certified_zero_late(1.0)  # type: ignore[arg-type]
