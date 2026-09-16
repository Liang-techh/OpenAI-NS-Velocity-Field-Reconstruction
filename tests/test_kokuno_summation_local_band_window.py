from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_summation_local_band_window import (
    FULL_RECONSTRUCTION,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    active_band_implication,
    dyadic_scale,
    every_active_band_is_in_window,
    local_band_window,
    verify_dyadic_indices,
)


def test_exact_boundary_case_realizes_source_span_four():
    window = local_band_window(Fraction(1, 64))
    assert window.lower_scale == Fraction(1, 256)
    assert window.upper_scale == Fraction(1, 16)
    assert window.indices == (4, 5, 6, 7, 8)
    assert window.scales == (
        Fraction(1, 16),
        Fraction(1, 32),
        Fraction(1, 64),
        Fraction(1, 128),
        Fraction(1, 256),
    )
    assert window.index_span == 4
    assert verify_dyadic_indices(window)


def test_non_dyadic_local_q_still_has_finite_exact_window():
    window = local_band_window(Fraction(3, 128))
    assert window.indices == (4, 5, 6, 7)
    assert window.index_span == 3
    assert verify_dyadic_indices(window)


def test_active_band_hypotheses_imply_conservative_window():
    q_local = Fraction(1, 64)
    samples = (
        (q_local, Fraction(1, 32)),
        (q_local, Fraction(1, 64)),
        (q_local, Fraction(1, 128)),
        (Fraction(3, 128), Fraction(1, 64)),
        (Fraction(9, 1024), Fraction(1, 128)),
    )
    assert every_active_band_is_in_window(q_local, samples)


def test_open_neighborhood_boundary_is_not_silently_admitted():
    q_local = Fraction(1, 64)
    assert not active_band_implication(q_local, q_local / 2, q_local / 2)
    assert not active_band_implication(q_local, 2 * q_local, q_local)


def test_wrong_band_scale_fails_instead_of_being_projected_into_window():
    q_local = Fraction(1, 64)
    q_point = q_local
    assert not active_band_implication(q_local, q_point, Fraction(1, 8))


def test_exact_dyadic_scale_supports_negative_indices():
    assert dyadic_scale(-3) == 8
    window = local_band_window(8)
    assert verify_dyadic_indices(window)
    assert window.index_span <= 4


def test_invalid_or_inexact_inputs_fail_closed():
    with pytest.raises(TypeError):
        local_band_window(0.125)
    with pytest.raises(TypeError):
        dyadic_scale(1.0)
    with pytest.raises(ValueError):
        local_band_window(0)
    with pytest.raises(ValueError):
        active_band_implication(1, 0, 1)
    with pytest.raises(RuntimeError):
        local_band_window(Fraction(1, 2**100), max_steps=10)


def test_status_flags_remain_false():
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
