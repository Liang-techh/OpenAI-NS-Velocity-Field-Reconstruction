from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_iteration_ledger import (
    PINNED_FORMAL_COMMIT,
    PINNED_KAPPA,
    Section9ResidualImprovementCertificate,
    certify_residual_improvement,
    gain,
    input_sigma,
    residual_mean,
    residual_minimum,
    residual_wave,
    sigma,
)


def test_pinned_schedule_is_exact_fraction_arithmetic() -> None:
    assert sigma(0) == Fraction(1, 5)
    assert sigma(1) == Fraction(3, 10)
    assert sigma(17) == Fraction(19, 10)
    assert input_sigma(1) == Fraction(1, 5)
    assert input_sigma(8) == sigma(7)
    assert gain(Fraction(3, 7), 4) == Fraction(6, 35)


def test_one_cycle_certifies_exact_native_and_physical_improvement() -> None:
    cert = certify_residual_improvement(7, Fraction(3, 11))

    assert cert.sigma_after - cert.sigma_before == Fraction(1, 10)
    assert cert.wave_after - cert.wave_before == Fraction(1, 10)
    assert cert.mean_after - cert.mean_before == Fraction(1, 10)
    assert cert.minimum_before == cert.wave_before
    assert cert.minimum_after == cert.wave_after
    assert cert.physical_increment == Fraction(3, 110)
    assert cert.physical_minimum_after - cert.physical_minimum_before == Fraction(3, 110)
    assert cert.kappa == Fraction(1, 100000)


def test_residual_minimum_is_wave_row_for_multiple_stages() -> None:
    for stage in (0, 1, 2, 19, 101):
        assert residual_minimum(stage) == residual_wave(stage)
        assert residual_wave(stage) < residual_mean(stage)


def test_truth_boundary_stays_finite_and_non_materialized() -> None:
    cert = certify_residual_improvement(0, 1)

    assert cert.finite_step_only is True
    assert cert.actual_stage_field_consumed is False
    assert cert.actual_residual_evaluated is False
    assert cert.infinite_correction_sequence_certified is False
    assert cert.eq_9_21_summed_field_certified is False
    assert cert.paper_exact_velocity_available is False


def test_float_and_unpinned_kappa_are_rejected() -> None:
    with pytest.raises(TypeError, match="exact integer/Fraction"):
        certify_residual_improvement(0, 0.1)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="exact integer/Fraction"):
        certify_residual_improvement(0, 1, kappa=1e-5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="1/100000"):
        certify_residual_improvement(0, 1, kappa=Fraction(2, 100000))


def test_invalid_stage_and_nonpositive_h_are_rejected() -> None:
    with pytest.raises(ValueError, match="nonnegative"):
        sigma(-1)
    with pytest.raises(TypeError, match="nonnegative integer"):
        sigma(True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="at least 1"):
        input_sigma(0)
    with pytest.raises(ValueError, match="strictly positive"):
        certify_residual_improvement(0, 0)


def test_pinned_revision_cross_wire_is_rejected() -> None:
    cert = Section9ResidualImprovementCertificate(stage=3, h=Fraction(2, 5))
    with pytest.raises(ValueError, match="formal_commit"):
        replace(cert, formal_commit=PINNED_FORMAL_COMMIT[:-1] + "0")


def test_pinned_kappa_constant_is_exact() -> None:
    assert PINNED_KAPPA == Fraction(1, 100000)
