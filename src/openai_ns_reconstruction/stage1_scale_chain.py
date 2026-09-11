"""Fail-closed bridge from the actual schedule to the Stage-1 scale chain.

The landed Stage-1 pieces construct an actual SchedulePressure analytic
neighborhood and, separately, turn coefficient-space norm ledgers into
AxisResolvent/remainder/Lambda-C bounds.  This module connects those pieces
without inventing fixed-point data.

The first version deliberately pushed the single common eleven-field bound
through every coefficient.  That reproduced the simplest pinned Lean existence
proof, but it also let the large gradient bound contaminate ``||chi||`` even
though ``AxisResolvent`` depends on chi specifically.  The pinned
``boundedAxisElement_norm`` theorem applies to one field and one value bound at
a time, so the actual-schedule path now uses the theorem-faithful componentwise
ledger while retaining the same common radius and a common family bound given
by the maximum of the individual norms.

Every numeric representation failure remains fail-closed.  Such an obstruction
is about this conservative certificate chain only; it is **not** a lower bound
on the true resolvent norm and is not a claim that the theorem's fixed point
fails to exist.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_FLOOR, localcontext
import math
import sys
from typing import Optional

from .axis_componentwise_input_bounds import (
    ComponentwiseAnalyticInputNormCertificate,
    componentwise_analytic_input_norm_certificate,
)
from .axis_remainder_bounds import NaturalRemainderBoundCertificate, natural_remainder_bound_certificate
from .natural_scale_selection import NaturalScaleSelection, select_natural_scale
from .outgoing_tail import TailData
from .schedule_analytic_neighborhood import (
    ActualScheduleAnalyticInputCertificate,
    certify_actual_schedule_analytic_inputs,
)

_RADIUS_LOSS_HALF = 12.0
_RESOLVENT_POWER_CONSTANT = 2560.0
_PROBE_PRECISION = 100
_DEFAULT_MAX_PROBE_TERMS = 256


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _mul_upper(left: float, right: float, name: str) -> float:
    left = _finite_positive(left, f"{name} left factor")
    right = _finite_positive(right, f"{name} right factor")
    value = left * right
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} overflowed binary64")
    return math.nextafter(value, math.inf)


def _is_finite_representation_error(exc: BaseException) -> bool:
    """Recognize downstream fail-closed errors caused by non-finite float bounds."""

    if isinstance(exc, (ArithmeticError, OverflowError)):
        return True
    if isinstance(exc, ValueError):
        text = str(exc).lower()
        return "must be finite" in text or "overflow" in text
    return False


@dataclass(frozen=True)
class Binary64MajorantObstruction:
    """One positive factorial-majorant term already exceeds binary64 range."""

    term_index: int
    majorant_parameter_upper: float
    term_lower: Decimal
    binary64_max: Decimal

    def __post_init__(self) -> None:
        if self.term_index < 0:
            raise ValueError("term_index must be nonnegative")
        _finite_positive(self.majorant_parameter_upper, "majorant_parameter_upper")
        if not self.term_lower.is_finite() or self.term_lower <= 0:
            raise ValueError("term_lower must be finite and positive")
        if self.term_lower <= self.binary64_max:
            raise ValueError("obstruction term must exceed the largest binary64 value")

    @property
    def decimal_order_lower(self) -> int:
        """Integer n such that the witnessed term is at least 10**n."""

        return int(self.term_lower.adjusted())


def first_binary64_majorant_obstruction(
    majorant_parameter_upper: float,
    *,
    max_terms: int = _DEFAULT_MAX_PROBE_TERMS,
) -> Optional[Binary64MajorantObstruction]:
    """Find a positive-term witness that the conservative series exceeds float.

    Decimal arithmetic is rounded *downward*.  Consequently each recurrence
    value is a lower bound for the corresponding positive term at the supplied
    binary64 ``K``.  If that lower bound exceeds ``sys.float_info.max``, the sum
    at this conservative ``K`` certainly cannot be converted to binary64.

    Returning ``None`` means only that no witness was found in this finite
    prefix; it does not certify that the complete series fits in binary64.
    """

    K = _finite_positive(majorant_parameter_upper, "majorant_parameter_upper")
    if not isinstance(max_terms, int) or max_terms < 1:
        raise ValueError("max_terms must be a positive integer")

    binary64_max = Decimal.from_float(sys.float_info.max)
    with localcontext() as ctx:
        ctx.prec = _PROBE_PRECISION
        ctx.rounding = ROUND_FLOOR
        kd = Decimal.from_float(K)
        term = Decimal(1)
        for k in range(1, max_terms + 1):
            term = (term * kd) / Decimal(k * (k + 1))
            if term > binary64_max:
                return Binary64MajorantObstruction(
                    term_index=k,
                    majorant_parameter_upper=K,
                    term_lower=+term,
                    binary64_max=binary64_max,
                )
    return None


@dataclass(frozen=True)
class Binary64PropagationObstruction:
    """A later certified bound cannot be represented by current float machinery."""

    stage: str
    detail: str

    def __post_init__(self) -> None:
        if self.stage not in {
            "resolvent-conversion",
            "remainder-propagation",
            "scale-selection",
        }:
            raise ValueError("unknown Stage-1 propagation obstruction")
        if not self.detail:
            raise ValueError("propagation obstruction detail must be nonempty")


@dataclass(frozen=True)
class ActualScheduleScaleChainDiagnostic:
    """Actual-schedule bridge through the first representable theorem layer."""

    schedule_inputs: ActualScheduleAnalyticInputCertificate
    coefficient_family_norm_upper: float
    resolvent_majorant_parameter_upper: float
    obstruction: Optional[Binary64MajorantObstruction]
    propagation_obstruction: Optional[Binary64PropagationObstruction]
    analytic_norms: Optional[ComponentwiseAnalyticInputNormCertificate]
    remainder: Optional[NaturalRemainderBoundCertificate]
    scale: Optional[NaturalScaleSelection]

    @property
    def reached_remainder_chain(self) -> bool:
        return self.remainder is not None

    @property
    def reached_scale_selection(self) -> bool:
        return self.scale is not None


def _stopped(
    *,
    schedule: ActualScheduleAnalyticInputCertificate,
    family_upper: float,
    K_upper: float,
    obstruction: Optional[Binary64MajorantObstruction] = None,
    propagation: Optional[Binary64PropagationObstruction] = None,
    analytic: Optional[ComponentwiseAnalyticInputNormCertificate] = None,
    remainder: Optional[NaturalRemainderBoundCertificate] = None,
) -> ActualScheduleScaleChainDiagnostic:
    return ActualScheduleScaleChainDiagnostic(
        schedule_inputs=schedule,
        coefficient_family_norm_upper=family_upper,
        resolvent_majorant_parameter_upper=K_upper,
        obstruction=obstruction,
        propagation_obstruction=propagation,
        analytic_norms=analytic,
        remainder=remainder,
        scale=None,
    )


def diagnose_actual_schedule_scale_chain(
    data: TailData,
    j: float,
    *,
    max_probe_terms: int = _DEFAULT_MAX_PROBE_TERMS,
) -> ActualScheduleScaleChainDiagnostic:
    """Connect actual schedule inputs to the landed norm/remainder/scale chain.

    The caller supplies only the theorem's actual ``TailData`` and ``j``.  No
    free ``rho``, per-field value bounds, resolvent norm, remainder constants,
    Lambda, or C are accepted here.

    Each fixed field now receives its own ``12*B_k`` Cauchy norm certificate.
    The common ``CoefficientFamily.bound`` is their maximum, while the pinned
    natural-resolvent majorant uses only the chi-specific norm.  Numeric
    overflow at any later layer is returned as a structured obstruction rather
    than being hidden or replaced by smaller sampled values.
    """

    schedule = certify_actual_schedule_analytic_inputs(data, j)
    neighborhood = schedule.neighborhood

    family_upper = _mul_upper(
        max(neighborhood.field_bounds.values()),
        _RADIUS_LOSS_HALF,
        "componentwise coefficient-family max norm bound",
    )
    chi_upper = _mul_upper(
        neighborhood.field_bounds["chi"],
        _RADIUS_LOSS_HALF,
        "chi coefficient norm bound",
    )
    K_upper = _mul_upper(
        _RESOLVENT_POWER_CONSTANT,
        chi_upper,
        "natural-resolvent majorant parameter",
    )

    obstruction = first_binary64_majorant_obstruction(
        K_upper,
        max_terms=max_probe_terms,
    )
    if obstruction is not None:
        return _stopped(
            schedule=schedule,
            family_upper=family_upper,
            K_upper=K_upper,
            obstruction=obstruction,
        )

    analytic = componentwise_analytic_input_norm_certificate(
        neighborhood_radius=neighborhood.radius,
        field_value_sup_upper=neighborhood.field_bounds,
    )
    if analytic.resolvent_majorant.majorant_parameter_upper != K_upper:
        raise ArithmeticError("componentwise chi ledger disagrees with Stage-1 K bound")

    try:
        operators = analytic.operator_norm_bounds()
    except Exception as exc:
        if not _is_finite_representation_error(exc):
            raise
        return _stopped(
            schedule=schedule,
            family_upper=family_upper,
            K_upper=K_upper,
            propagation=Binary64PropagationObstruction(
                "resolvent-conversion",
                str(exc),
            ),
            analytic=analytic,
        )

    axis_data = analytic.axis_data_norm_bounds(h=data.h)
    try:
        remainder = natural_remainder_bound_certificate(
            operators=operators,
            data=axis_data,
            amplitude_norm_upper=analytic.amplitude_norm_upper,
        )
    except Exception as exc:
        if not _is_finite_representation_error(exc):
            raise
        return _stopped(
            schedule=schedule,
            family_upper=family_upper,
            K_upper=K_upper,
            propagation=Binary64PropagationObstruction(
                "remainder-propagation",
                str(exc),
            ),
            analytic=analytic,
        )

    try:
        scale = select_natural_scale(
            remainder_bound=remainder.remainder_bound_upper,
            remainder_lipschitz=remainder.remainder_lipschitz_upper,
            phase_real_part_sup=neighborhood.axis_phase_real_part_sup_upper,
        )
    except Exception as exc:
        if not _is_finite_representation_error(exc):
            raise
        return _stopped(
            schedule=schedule,
            family_upper=family_upper,
            K_upper=K_upper,
            propagation=Binary64PropagationObstruction(
                "scale-selection",
                str(exc),
            ),
            analytic=analytic,
            remainder=remainder,
        )

    return ActualScheduleScaleChainDiagnostic(
        schedule_inputs=schedule,
        coefficient_family_norm_upper=family_upper,
        resolvent_majorant_parameter_upper=K_upper,
        obstruction=None,
        propagation_obstruction=None,
        analytic_norms=analytic,
        remainder=remainder,
        scale=scale,
    )
