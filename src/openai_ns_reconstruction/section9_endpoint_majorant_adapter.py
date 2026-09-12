"""Analytic adapter from a Section 9 Eq. (9.18) envelope to an endpoint majorant.

Section 10 needs an integrable physical-time bound for ``partial_t D^n R`` on
its fixed late-time support.  The generic endpoint bridge intentionally accepts
such a bound only as theorem input.  This module removes one avoidable free
parameter: when a theorem source supplies a *uniform* Eq. (9.18) envelope on
the genuine Section 9 residual norm, the coefficient of a bounded
``(1-t)^0`` endpoint majorant is derived here rather than supplied independently.

On the official Section 10 late plateau ``t >= 3/4`` and support
``|z| <= 1/4``, Eq. (4.1) gives

    0 < q <= (1-t) + z^2 <= 5/16.

For beta > 0,

    q^beta (1 + |log q|)^P

has a finite, explicitly computable supremum on ``0 < q <= 5/16``.  If
``beta = 0`` the adapter only accepts ``P = 0``.  Negative beta, or
``beta = 0, P > 0``, is rejected because this route does not give a uniform
endpoint bound.

The input remains theorem-facing: this file does not construct the Section 9
correction sequence, prove Eq. (9.18), or certify that the residual provider is
the manuscript one.  It only converts a source-bound uniform envelope into the
already-landed Section 9 -> Section 10 endpoint-majorant interface.  Therefore
all paper-exact truth flags remain false.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral, Rational, Real

from .section9_endpoint_bridge import (
    Section9EndpointBridgeRecord,
    Section9UniformEndpointDerivativeWitness,
    admit_section9_uniform_endpoint_witness,
)
from .section9_physical_window import (
    LATE_START_EXACT,
    SECTION10_ENDPOINT_EXACT,
    SUPPORT_Z_SQUARED_EXACT,
)
from .section9_residual_decay import Scalar, section9_residual_gain
from .section9_stage_certificate import CertifiedBoundDatum


OFFICIAL_SECTION10_Q_MAX = (
    SECTION10_ENDPOINT_EXACT - LATE_START_EXACT + SUPPORT_Z_SQUARED_EXACT
)


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _fraction(value: object, name: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real number")
    if isinstance(value, Fraction):
        out = value
    elif isinstance(value, Integral):
        out = Fraction(int(value), 1)
    elif isinstance(value, Rational):
        out = Fraction(value)
    elif isinstance(value, Real):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f"{name} must be finite")
        out = Fraction(str(numeric))
    else:
        raise TypeError(f"{name} must be a real rationalizable value")
    return out


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _power_log_supremum(beta: Fraction, log_power: int) -> tuple[float, float]:
    """Return ``(supremum, maximizing_s)`` for the official q interval.

    With ``s=-log(q)``, the factor is ``exp(-beta*s) * (1+s)^P`` on
    ``s >= -log(5/16)``.  For positive beta the log derivative is
    ``-beta + P/(1+s)``, so the maximum is attained at the left endpoint or
    the unique stationary point.  The beta=0, P=0 constant case is also
    admitted.
    """

    log_power = _natural(log_power, "log_power")
    if beta < 0 or (beta == 0 and log_power > 0):
        raise ValueError(
            "Eq. (9.18) leading exponent must be positive, except beta=0 is "
            "allowed only when log_power=0"
        )

    s0 = -math.log(float(OFFICIAL_SECTION10_Q_MAX))
    if beta == 0:
        return 1.0, s0

    beta_float = float(beta)
    stationary = (log_power / beta_float) - 1.0 if log_power else -math.inf
    s_star = max(s0, stationary)
    log_sup = -beta_float * s_star + log_power * math.log1p(s_star)
    sup = math.exp(log_sup)
    if not math.isfinite(sup) or sup <= 0.0:
        raise ArithmeticError("failed to derive a finite positive power-log supremum")
    return sup, s_star


@dataclass(frozen=True)
class Section9UniformResidualEnvelopeWitness:
    """Theorem input for a uniform Eq. (9.18) envelope on the endpoint domain.

    ``leading_constant`` is ``C_{j,m}`` and ``flat_remainder_uniform_bound`` is
    a uniform upper bound for ``E_{j,m}`` on the stated source/domain.  Both
    must come from the same evidence class.  The three boolean fields are
    explicit theorem assertions: this adapter cannot infer them from sampled
    values or from provenance strings.

    ``valid_q_upper`` is the strict upper end of the q-domain where the theorem
    envelope applies.  It must exceed the analytically derived official
    Section 10 maximum ``5/16``.
    """

    endpoint_derivative_degree: int
    stage: int
    spatial_window: int
    h: Scalar
    derivative_loss: Scalar
    log_power: int
    leading_constant: CertifiedBoundDatum
    flat_remainder_uniform_bound: CertifiedBoundDatum
    valid_q_upper: Fraction | float | int
    source_id: str
    source_revision: str
    residual_envelope_uniform_certified: bool
    flat_remainder_uniform_certified: bool
    official_support_covered_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "endpoint_derivative_degree",
            _natural(self.endpoint_derivative_degree, "endpoint_derivative_degree"),
        )
        object.__setattr__(self, "stage", _natural(self.stage, "stage"))
        object.__setattr__(
            self, "spatial_window", _natural(self.spatial_window, "spatial_window")
        )
        object.__setattr__(self, "log_power", _natural(self.log_power, "log_power"))
        if not isinstance(self.leading_constant, CertifiedBoundDatum):
            raise TypeError("leading_constant must be a CertifiedBoundDatum")
        if not isinstance(self.flat_remainder_uniform_bound, CertifiedBoundDatum):
            raise TypeError(
                "flat_remainder_uniform_bound must be a CertifiedBoundDatum"
            )
        if self.leading_constant.upper_bound <= 0.0:
            raise ValueError("Eq. (9.18) leading constant must be positive")
        if self.leading_constant.kind != self.flat_remainder_uniform_bound.kind:
            raise ValueError(
                "leading and flat-remainder bounds must use the same evidence kind"
            )

        valid_q_upper = _fraction(self.valid_q_upper, "valid_q_upper")
        if not OFFICIAL_SECTION10_Q_MAX < valid_q_upper <= 1:
            raise ValueError(
                "valid_q_upper must satisfy 5/16 < valid_q_upper <= 1 so the "
                "uniform Eq. (9.18) envelope covers the official Section 10 support"
            )
        object.__setattr__(self, "valid_q_upper", valid_q_upper)
        object.__setattr__(
            self, "source_id", _nonempty_text(self.source_id, "source_id")
        )
        object.__setattr__(
            self,
            "source_revision",
            _nonempty_text(self.source_revision, "source_revision"),
        )
        for name in (
            "residual_envelope_uniform_certified",
            "flat_remainder_uniform_certified",
            "official_support_covered_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def section9_derivative_order(self) -> int:
        return self.endpoint_derivative_degree + 1

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision


@dataclass(frozen=True)
class Section9DerivedEndpointMajorantRecord:
    """Machine-derived bounded endpoint majorant from one theorem envelope."""

    source_id: str
    source_revision: str
    endpoint_derivative_degree: int
    section9_derivative_order: int
    residual_exponent: Fraction
    official_q_max: Fraction
    power_log_supremum: float
    maximizing_log_coordinate: float
    derived_uniform_coefficient: float
    bridge: Section9EndpointBridgeRecord
    status: str = "formal-structure"
    actual_section9_sequence_verified: bool = False
    source_majorant_derived_from_actual_residual_verified: bool = False
    endpoint_limit_constructed: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def formal_adapter_ready(self) -> bool:
        return (
            self.section9_derivative_order == self.endpoint_derivative_degree + 1
            and self.official_q_max == OFFICIAL_SECTION10_Q_MAX
            and self.power_log_supremum > 0.0
            and math.isfinite(self.power_log_supremum)
            and self.derived_uniform_coefficient >= 0.0
            and math.isfinite(self.derived_uniform_coefficient)
            and self.bridge.formal_bridge_ready
            and self.bridge.endpoint_derivative_degree
            == self.endpoint_derivative_degree
            and self.bridge.section9_derivative_order
            == self.section9_derivative_order
            and self.bridge.majorant.singularity_exponent == 0.0
            and self.bridge.majorant.valid_from == float(LATE_START_EXACT)
            and not self.actual_section9_sequence_verified
            and not self.source_majorant_derived_from_actual_residual_verified
            and not self.endpoint_limit_constructed
            and not self.paper_exact_velocity_available
        )


def derive_section9_endpoint_majorant_from_uniform_envelope(
    witness: Section9UniformResidualEnvelopeWitness,
) -> Section9DerivedEndpointMajorantRecord:
    """Derive an ``alpha=0`` endpoint majorant from a uniform Eq. (9.18) envelope.

    The Section 9 residual norm of order ``m=n+1`` controls
    ``partial_t D^n R``.  The official physical support gives ``q<=5/16``.
    This function takes the exact Section 9 exponent

        beta = h*sigma_j - K_m

    and analytically maximizes the power-log factor on that q interval.  The
    flat remainder is then added as an independently supplied *uniform*
    theorem bound.  No coefficient/exponent fitting or residual sampling is
    performed.
    """

    if not isinstance(witness, Section9UniformResidualEnvelopeWitness):
        raise TypeError("witness must be a Section9UniformResidualEnvelopeWitness")

    derivative_loss = _fraction(witness.derivative_loss, "derivative_loss")
    if derivative_loss < 0:
        raise ValueError("derivative_loss must be nonnegative")

    residual_exponent = section9_residual_gain(witness.stage, witness.h) - derivative_loss
    factor, s_star = _power_log_supremum(residual_exponent, witness.log_power)

    derived = (
        witness.leading_constant.upper_bound * factor
        + witness.flat_remainder_uniform_bound.upper_bound
    )
    if not math.isfinite(derived) or derived < 0.0:
        raise ArithmeticError("derived endpoint coefficient is not finite and nonnegative")

    provenance = (
        f"{witness.source_id}@{witness.source_revision}; "
        f"Eq. (9.18) uniform leading envelope: {witness.leading_constant.provenance}; "
        "uniform flat remainder: "
        f"{witness.flat_remainder_uniform_bound.provenance}; "
        "coefficient analytically derived on official Section 10 "
        "t>=3/4, |z|<=1/4 support with q<=5/16"
    )
    coefficient = CertifiedBoundDatum(
        upper_bound=derived,
        kind=witness.leading_constant.kind,
        provenance=provenance,
    )
    endpoint_witness = Section9UniformEndpointDerivativeWitness(
        endpoint_derivative_degree=witness.endpoint_derivative_degree,
        section9_derivative_order=witness.section9_derivative_order,
        stage=witness.stage,
        spatial_window=witness.spatial_window,
        coefficient_bound=coefficient,
        singularity_exponent=0.0,
        valid_from=float(LATE_START_EXACT),
    )
    bridge = admit_section9_uniform_endpoint_witness(endpoint_witness)
    result = Section9DerivedEndpointMajorantRecord(
        source_id=witness.source_id,
        source_revision=witness.source_revision,
        endpoint_derivative_degree=witness.endpoint_derivative_degree,
        section9_derivative_order=witness.section9_derivative_order,
        residual_exponent=residual_exponent,
        official_q_max=OFFICIAL_SECTION10_Q_MAX,
        power_log_supremum=factor,
        maximizing_log_coordinate=s_star,
        derived_uniform_coefficient=derived,
        bridge=bridge,
    )
    if not result.formal_adapter_ready:
        raise ArithmeticError("Section 9 endpoint-majorant adapter invariant failed")
    return result
