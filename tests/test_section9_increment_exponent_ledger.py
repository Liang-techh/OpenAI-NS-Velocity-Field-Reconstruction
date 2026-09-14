from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_increment_exponent_ledger import (
    PINNED_KAPPA,
    Section9IncrementExponentCertificate,
    coordinate_A,
    gain,
    mean_native,
    radial_native,
    wave_native,
    wave_pressure_native,
)


def test_exact_native_formulas_and_successor_gain():
    c = Section9IncrementExponentCertificate(stage=7, h=Fraction(1, 4))
    assert wave_native(PINNED_KAPPA, 7) == Fraction(7, 10) + Fraction(3, 5) - PINNED_KAPPA
    assert wave_pressure_native(PINNED_KAPPA, 7) == Fraction(7, 10) + Fraction(11, 10) - PINNED_KAPPA
    assert mean_native(PINNED_KAPPA, 7) == Fraction(7, 10) + Fraction(11, 10) - 2 * PINNED_KAPPA
    assert radial_native(PINNED_KAPPA, 7) == mean_native(PINNED_KAPPA, 7) + 1
    assert set(c.native_successor_increments.values()) == {Fraction(1, 10)}
    assert c.common_gain == gain(Fraction(1, 4), 7) == Fraction(7, 40)
    assert c.common_gain_successor_increment == Fraction(1, 40)


def test_fixed_offsets_cancel_exactly_and_leave_positive_slack():
    h = Fraction(1, 4)
    c = Section9IncrementExponentCertificate(stage=7, h=h)
    assert coordinate_A(h) == Fraction(3, 4)
    assert c.fixed_offsets == {
        "wavePotential": h,
        "meanStream": 0,
        "directAngular": 0,
        "wavePressure": Fraction(3, 2),
        "meanPressure": 0,
    }
    assert c.physical_component_exponents == {
        "wavePotential": h * wave_native(PINNED_KAPPA, 7),
        "meanStream": h * mean_native(PINNED_KAPPA, 7),
        "directAngular": h * mean_native(PINNED_KAPPA, 7),
        "wavePressure": h * wave_pressure_native(PINNED_KAPPA, 7),
        "meanPressure": h * mean_native(PINNED_KAPPA, 7),
    }
    assert c.common_gain_lower_bound_holds
    assert c.component_slacks == {
        "wavePotential": h * (Fraction(3, 5) - PINNED_KAPPA),
        "meanStream": h * (Fraction(11, 10) - 2 * PINNED_KAPPA),
        "directAngular": h * (Fraction(11, 10) - 2 * PINNED_KAPPA),
        "wavePressure": h * (Fraction(11, 10) - PINNED_KAPPA),
        "meanPressure": h * (Fraction(11, 10) - 2 * PINNED_KAPPA),
    }
    assert min(c.component_slacks.values()) > 0


@pytest.mark.parametrize("bad_h", [0, Fraction(1, 2), Fraction(3, 4)])
def test_paper_h_regime_is_fail_closed(bad_h):
    with pytest.raises(ValueError, match="0 < h < 1/2"):
        Section9IncrementExponentCertificate(stage=1, h=bad_h)


@pytest.mark.parametrize("bad_stage", [0, -1])
def test_positive_increment_index_required(bad_stage):
    with pytest.raises(ValueError, match="at least 1"):
        Section9IncrementExponentCertificate(stage=bad_stage, h=Fraction(1, 4))


@pytest.mark.parametrize("bad_exact", [0.25, Decimal("0.25")])
def test_approximate_theorem_inputs_are_rejected(bad_exact):
    with pytest.raises(TypeError, match="exact integer/Fraction"):
        Section9IncrementExponentCertificate(stage=1, h=bad_exact)


def test_kappa_and_provenance_are_pinned():
    with pytest.raises(ValueError, match="1/100000"):
        Section9IncrementExponentCertificate(
            stage=1, h=Fraction(1, 4), kappa=Fraction(1, 99999)
        )
    with pytest.raises(ValueError, match="formal_commit"):
        Section9IncrementExponentCertificate(
            stage=1,
            h=Fraction(1, 4),
            formal_commit="not-the-pinned-revision",
        )


def test_truth_boundary_remains_finite_and_nonmaterialized():
    c = Section9IncrementExponentCertificate(stage=3, h=Fraction(1, 5))
    assert c.exact_increment_ledger_machine_checked
    assert c.finite_stage_only
    assert not c.actual_stage_field_consumed
    assert not c.support_estimates_consumed
    assert not c.infinite_correction_sequence_certified
    assert not c.eq_9_21_summed_field_certified
    assert not c.paper_exact_velocity_available
