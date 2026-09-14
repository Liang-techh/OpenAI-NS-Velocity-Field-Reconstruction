from fractions import Fraction

import pytest

from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.stage1_actual_chi_cutoff import (
    CHI_EXACT_LOWER,
    CHI_REQUIRED_LOWER,
    CHI_STRICT_SLACK,
    SIGMA_SQUARED_MARGIN_RATIO,
    ActualScheduleChiCutoffCertificate,
    actual_schedule_chi_cutoff_certificate,
)


def _data(*, P: float = 2.0, h: float = 0.01) -> TailData:
    return TailData(
        OutgoingCoreParameters(P=P, m=1.0, lam=0.05, wait=30.0),
        h=h,
    )


def test_actual_schedule_margin_closes_exact_chi_threshold_algebra() -> None:
    certificate = actual_schedule_chi_cutoff_certificate(_data(), 0.05)

    assert certificate.margin_witness.margin > 0.0
    assert certificate.low_Z_delta == certificate.margin_witness.parameters.j / 10.0
    assert certificate.sigma == certificate.margin_witness.cutoff.sigma
    assert certificate.sigma_squared_margin_ratio == Fraction(1, 400)
    assert certificate.chi_lower_exact == Fraction(400, 401)
    assert certificate.chi_required_lower == Fraction(99, 100)
    assert certificate.strict_slack == Fraction(301, 40100)
    assert CHI_STRICT_SLACK == Fraction(301, 40100)
    assert certificate.chi_required_lower < certificate.chi_lower_exact
    assert certificate.implication == "|Z(eta)| <= j/10 -> 99/100 < chi(h,j,sigma,eta)"
    assert certificate.paper_exact is False
    assert certificate.full_reconstruction is False


def test_downstream_chi_admission_requires_exact_rational_evidence() -> None:
    certificate = actual_schedule_chi_cutoff_certificate(_data(), 0.05)

    certificate.require_exact_downstream_bound(Fraction(99, 100))
    certificate.require_exact_downstream_bound(Fraction(400, 401))

    with pytest.raises(TypeError, match="exact Fraction"):
        certificate.require_exact_downstream_bound(0.99)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="below"):
        certificate.require_exact_downstream_bound(Fraction(98, 100))
    with pytest.raises(ValueError, match="exceeds"):
        certificate.require_exact_downstream_bound(Fraction(401, 401))


def test_certificate_rejects_drifted_exact_algebra() -> None:
    base = actual_schedule_chi_cutoff_certificate(_data(), 0.05)

    with pytest.raises(ValueError, match="1/400"):
        ActualScheduleChiCutoffCertificate(
            margin_witness=base.margin_witness,
            sigma_squared_margin_ratio=Fraction(1, 399),
            chi_lower_exact=CHI_EXACT_LOWER,
            chi_required_lower=CHI_REQUIRED_LOWER,
        )

    with pytest.raises(ValueError, match="400/401"):
        ActualScheduleChiCutoffCertificate(
            margin_witness=base.margin_witness,
            sigma_squared_margin_ratio=SIGMA_SQUARED_MARGIN_RATIO,
            chi_lower_exact=Fraction(399, 400),
            chi_required_lower=CHI_REQUIRED_LOWER,
        )


def test_actual_schedule_constructor_fails_closed_outside_margin_hypotheses() -> None:
    with pytest.raises(ValueError, match="P >= 2"):
        actual_schedule_chi_cutoff_certificate(_data(P=1.5), 0.05)

    with pytest.raises(ValueError, match="1/100"):
        actual_schedule_chi_cutoff_certificate(_data(h=0.02), 0.05)

    with pytest.raises(ValueError, match="1/20"):
        actual_schedule_chi_cutoff_certificate(_data(), 0.051)

    with pytest.raises(TypeError, match="TailData"):
        actual_schedule_chi_cutoff_certificate(object(), 0.05)  # type: ignore[arg-type]
