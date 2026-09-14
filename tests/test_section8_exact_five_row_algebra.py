from fractions import Fraction
from decimal import Decimal

import pytest

from openai_ns_reconstruction.section8_exact_five_row_algebra import (
    ExactFiveRowAlgebraCertificate,
    PINNED_FORMAL_COMMIT,
    PINNED_FIVE_ROWS_SYMBOL,
    PINNED_NORMALIZE_DEBT_SYMBOL,
    PINNED_PHYSICAL_FIVE_ROWS_SYMBOL,
)


def _certificate() -> ExactFiveRowAlgebraCertificate:
    return ExactFiveRowAlgebraCertificate(
        lam=Fraction(1, 2),
        C=2,
        ell=3,
        U=2,
        debt=(8, 108, 72),
    )


def test_exact_five_row_targets_and_zero_auxiliary_averages():
    certificate = _certificate()

    assert certificate.normalized_debt == (
        Fraction(2),
        Fraction(1),
        Fraction(2),
    )
    assert certificate.angular_powers == (
        Fraction(2),
        Fraction(-3),
        Fraction(-1),
    )
    assert certificate.axial_powers == (Fraction(1), Fraction(0))
    assert certificate.angular_moment_targets == (
        Fraction(0),
        Fraction(-1, 2),
        Fraction(1),
    )
    assert certificate.axial_moment_targets == (Fraction(0), Fraction(-1, 2))
    assert certificate.zero_auxiliary_average_identity_holds is True
    assert certificate.certified_physical_rows == (
        Fraction(0),
        Fraction(0),
        Fraction(-8),
        Fraction(-108),
        Fraction(-72),
    )
    assert certificate.certified_physical_rows == certificate.physical_row_targets
    assert certificate.five_row_identity_holds is True


def test_exact_row_map_detects_a_nonzero_defect_without_tolerance():
    certificate = _certificate()
    altered_angular = list(certificate.angular_moment_targets)
    altered_angular[1] += Fraction(1, 10**30)

    rows = certificate.physical_rows_from_normalized_moments(
        altered_angular,
        certificate.axial_moment_targets,
    )

    assert rows != certificate.physical_row_targets
    assert rows[2] - certificate.physical_row_targets[2] == Fraction(16, 10**30)


def test_exact_certificate_rejects_approximate_or_invalid_inputs():
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        ExactFiveRowAlgebraCertificate(0.5, 2, 3, 2, (8, 108, 72))
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        ExactFiveRowAlgebraCertificate(Fraction(1, 2), Decimal("2"), 3, 2, (8, 108, 72))
    with pytest.raises(ValueError, match="lam must be positive"):
        ExactFiveRowAlgebraCertificate(0, 2, 3, 2, (8, 108, 72))
    with pytest.raises(ValueError, match="C must be nonzero"):
        ExactFiveRowAlgebraCertificate(Fraction(1, 2), 0, 3, 2, (8, 108, 72))
    with pytest.raises(ValueError, match="ell must be positive"):
        ExactFiveRowAlgebraCertificate(Fraction(1, 2), 2, 0, 2, (8, 108, 72))
    with pytest.raises(ValueError, match="U must be nonzero"):
        ExactFiveRowAlgebraCertificate(Fraction(1, 2), 2, 3, 0, (8, 108, 72))


def test_truth_boundary_and_pinned_provenance_are_explicit():
    certificate = _certificate()

    assert PINNED_FORMAL_COMMIT == "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
    assert PINNED_FIVE_ROWS_SYMBOL == "NavierStokes.FiveRowRank.FiveRows"
    assert PINNED_NORMALIZE_DEBT_SYMBOL == "NavierStokes.MeanRankUpdate.normalizeDebt"
    assert PINNED_PHYSICAL_FIVE_ROWS_SYMBOL == "NavierStokes.MeanRankUpdate.physical_five_rows"
    assert certificate.status == "formal-structure"
    assert certificate.exact_five_row_algebra_machine_checked is True
    assert certificate.localized_moment_inverse_materialized is False
    assert certificate.actual_wave_defect_consumed is False
    assert certificate.actual_compact_correction_field_materialized is False
    assert certificate.paper_exact_velocity_available is False
