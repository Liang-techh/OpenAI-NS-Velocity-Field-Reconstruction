from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.natural_scale_selection_wide import (
    select_natural_scale_wide,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.stage1_natural_entrance_certificate import (
    CHI_LOG_SLOPE_MIN,
    ENTRANCE_SCALED_RADIUS_MAX,
    LOG_SLOPE_SCALED_RADIUS,
    LOG_SLOPE_STRICT_LOWER,
    PHI_STRICT_LOWER,
    NaturalEntranceScaleCertificate,
    actual_schedule_natural_entrance_scale_certificate,
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_schedule_closes_only_the_natural_entrance_scalar_gate() -> None:
    cert = actual_schedule_natural_entrance_scale_certificate(_schedule_data(), 0.05)

    assert cert.Lambda >= cert.stability_threshold
    assert cert.fixed_point_error_radius_safe <= Decimal(1)
    assert cert.remainder_bound_upper.is_finite()
    assert cert.paper_exact is False
    assert cert.full_reconstruction is False

    # These are exact formal theorem constants, not sampled profile claims.
    assert ENTRANCE_SCALED_RADIUS_MAX == Fraction(41, 10)
    assert PHI_STRICT_LOWER == Fraction(1, 8)
    assert LOG_SLOPE_SCALED_RADIUS == Fraction(4, 1)
    assert CHI_LOG_SLOPE_MIN == Fraction(99, 100)
    assert LOG_SLOPE_STRICT_LOWER == Fraction(23, 10)

    # The safe lower radius is admissible. No coefficient value is
    # manufactured by this regression.
    cert.require_backend_norm_error(cert.fixed_point_error_radius_safe)


def test_scale_certificate_rejects_nonpinned_stability_metadata() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(2),
        remainder_lipschitz=Decimal(3),
        phase_real_part_sup=Decimal("0.25"),
    )
    broken = type(scale)(
        remainder_bound=scale.remainder_bound,
        remainder_lipschitz=scale.remainder_lipschitz,
        phase_real_part_sup=scale.phase_real_part_sup,
        contraction=scale.contraction,
        stability=scale.stability + Decimal(1),
        Lambda=scale.Lambda + Decimal(1),
        C=scale.C,
    )

    with pytest.raises(ValueError, match="does not equal the pinned"):
        NaturalEntranceScaleCertificate.from_scale(broken)


def test_scale_certificate_rejects_lambda_below_stability_threshold() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(2),
        remainder_lipschitz=Decimal(0),
        phase_real_part_sup=Decimal(0),
    )
    broken = type(scale)(
        remainder_bound=scale.remainder_bound,
        remainder_lipschitz=scale.remainder_lipschitz,
        phase_real_part_sup=scale.phase_real_part_sup,
        contraction=scale.contraction,
        stability=scale.stability,
        Lambda=scale.stability - Decimal(1),
        C=scale.C,
    )

    with pytest.raises(ValueError, match="below its pinned stability"):
        NaturalEntranceScaleCertificate.from_scale(broken)


def test_backend_norm_error_gate_is_exact_decimal_and_fail_closed() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(2),
        remainder_lipschitz=Decimal(3),
        phase_real_part_sup=Decimal(0),
    )
    cert = NaturalEntranceScaleCertificate.from_scale(scale)

    cert.require_backend_norm_error(cert.fixed_point_error_radius_safe)

    with pytest.raises(ValueError, match=r"exceeds safe K/\(2\*Lambda\) radius"):
        cert.require_backend_norm_error(cert.fixed_point_error_radius_safe + Decimal(1))
    with pytest.raises(ValueError, match="finite nonnegative Decimal"):
        cert.require_backend_norm_error(0.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="finite nonnegative Decimal"):
        cert.require_backend_norm_error(Decimal("NaN"))


def test_log_slope_chi_gate_accepts_only_exact_theorem_bounds() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(1),
        remainder_lipschitz=Decimal(1),
        phase_real_part_sup=Decimal(0),
    )
    cert = NaturalEntranceScaleCertificate.from_scale(scale)

    cert.require_log_slope_chi_bounds(
        chi_lower=Fraction(99, 100),
        chi_upper=Fraction(1, 1),
    )

    with pytest.raises(ValueError, match="below 99/100"):
        cert.require_log_slope_chi_bounds(
            chi_lower=Fraction(989, 1000),
            chi_upper=Fraction(1, 1),
        )
    with pytest.raises(ValueError, match="exceeds 1"):
        cert.require_log_slope_chi_bounds(
            chi_lower=Fraction(99, 100),
            chi_upper=Fraction(1001, 1000),
        )
    with pytest.raises(TypeError, match="exact Fraction"):
        cert.require_log_slope_chi_bounds(
            chi_lower=0.99,  # type: ignore[arg-type]
            chi_upper=Fraction(1, 1),
        )
