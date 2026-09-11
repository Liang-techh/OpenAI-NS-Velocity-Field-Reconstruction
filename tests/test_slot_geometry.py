from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.charts import DyadicChart
from openai_ns_reconstruction.slow_labels import SlowLabel
from openai_ns_reconstruction.slot_geometry import (
    COVER_NORM_BOUND,
    PALETTE_CARD,
    ExplicitSlotSystem,
    covering_power,
    native_gap,
    palette_coordinates,
    palette_index,
    torus_sup_distance,
)


def _manual_matmul(a, b):
    return (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0],
         a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0],
         a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )


def _manual_cover_power(n):
    out = ((1, 0), (0, 1))
    for _ in range(n):
        out = _manual_matmul(((3, 1), (1, 5)), out)
    return out


def test_palette_matches_explicit_mod_9_mod_5_sign_construction():
    label = SlowLabel(ell=101, a=(-7, 8, 14), sigma=-1)
    assert palette_coordinates(label) == (2, (3, 3, 4), 0)
    idx = palette_index(label)
    assert 0 <= idx < PALETTE_CARD

    opposite = label.opposite()
    assert palette_index(opposite) == idx + 1

    # At equal level, an enlarged-box grid displacement <=4 cannot alias mod 5.
    shifted = SlowLabel(ell=101, a=(-3, 8, 14), sigma=-1)
    assert palette_index(shifted) != idx


def test_native_gap_and_cover_power_crosscheck_landed_charts():
    h = 0.005
    assert native_gap(h) == 5
    for n in range(6):
        assert covering_power(n) == _manual_cover_power(n)

    # This is the executable counterpart of SlotColoring.nativeIndex_gap for
    # representative landed chart levels, not an assumed equality in the code.
    D = native_gap(h)
    for ell in (20, 50, 100, 250):
        i = DyadicChart(ell, h).covering_index
        for other in range(ell - 4, ell + 5):
            j = DyadicChart(other, h).covering_index
            assert abs(i - j) <= D


def test_common_radius_has_exact_global_safety_margins():
    system = ExplicitSlotSystem(h=0.005)
    six_d = COVER_NORM_BOUND ** system.gap

    assert system.denominator == (PALETTE_CARD + 1) * six_d
    assert system.universal_center_gap == Fraction(1, system.denominator)
    assert system.rectangle_radius == system.standard_radius / 3

    # Independent algebraic checks of the two global inequalities used to
    # choose the conservative standard radius.
    assert 4 * system.standard_radius * six_d < 1
    assert 2 * system.standard_radius * (six_d + 1) < system.universal_center_gap
    assert math.isfinite(system.rectangle_radius_float)
    assert system.rectangle_radius_float > 0.0


def test_rational_centers_have_nonzero_forbidden_torus_gap_on_small_sample():
    system = ExplicitSlotSystem(h=0.005)
    colors = (0, 17, 2249)
    den = system.denominator

    # Independently recompute several forbidden center relations with the
    # manual integer matrix powers. Production does not enumerate these pairs.
    for n in range(system.gap + 1):
        matrix = _manual_cover_power(n)
        for i in colors:
            for j in colors:
                if n == 0 and i == j:
                    continue
                ci = (Fraction(i + 1, den), Fraction(0, 1))
                cj = (Fraction(j + 1, den), Fraction(0, 1))
                covered = (
                    matrix[0][0] * ci[0] + matrix[0][1] * ci[1],
                    matrix[1][0] * ci[0] + matrix[1][1] * ci[1],
                )
                assert torus_sup_distance(covered, cj) >= Fraction(1, den)


def test_pair_certificate_covers_sign_duplicate_same_level_and_cross_level():
    system = ExplicitSlotSystem(h=0.005)

    plus = SlowLabel(ell=100, a=(12, -3, 8), sigma=1)
    minus = plus.opposite()
    sign_cert = system.pair_certificate(plus, minus)
    assert sign_cert.certified
    assert sign_cert.covering_gap == 0
    assert sign_cert.lower_color != sign_cert.upper_color

    # Cross-level interaction needs no same-grid comparison: level mod 9
    # already separates all levels at distance <=4.
    left = SlowLabel(ell=100, a=(10_000, -20_000, 7), sigma=1)
    right = SlowLabel(ell=103, a=(-99_000, 88_000, -4), sigma=1)
    cross_cert = system.pair_certificate(left, right)
    assert cross_cert.certified
    assert cross_cert.covering_gap <= system.gap
    assert cross_cert.separation_perturbation_bound < cross_cert.universal_center_gap
    assert cross_cert.injection_perturbation_bound < 1


def test_interaction_certificate_fails_closed_outside_proved_discrete_hypotheses():
    system = ExplicitSlotSystem(h=0.005)
    base = SlowLabel(ell=100, a=(0, 0, 0), sigma=1)

    with pytest.raises(ValueError, match="distinct"):
        system.pair_certificate(base, base)

    with pytest.raises(ValueError, match="dyadic-level gap"):
        system.pair_certificate(base, SlowLabel(ell=105, a=(0, 0, 0), sigma=-1))

    # At equal level a displacement of five can alias the mod-5 grid colour;
    # physical overlap must rule this case out before the colouring theorem is used.
    with pytest.raises(ValueError, match="grid gaps"):
        system.pair_certificate(base, SlowLabel(ell=100, a=(5, 0, 0), sigma=1))

    with pytest.raises(ValueError, match="color"):
        system.center_for_color(PALETTE_CARD)
