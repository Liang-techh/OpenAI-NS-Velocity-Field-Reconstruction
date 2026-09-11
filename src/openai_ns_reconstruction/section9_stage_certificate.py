"""Fail-closed pointwise checks for the Section 9 stage estimates.

Lemma 9.8 gives, for every Cartesian space-time derivative order ``m``, the
stage-``j`` bounds

    |Z_j|_m + |Delta u_j|_m
      <= C_{j,m} q^(g_j-ell_m) (1 + |log q|)^P_{j,m},             (9.17)

and

    |R(u^[j], p^[j])|_m
      <= C_{j,m} q^(h sigma_j-K_m) (1 + |log q|)^P_{j,m}
         + E_{j,m}.                                               (9.18)

The exact exponents already live in :mod:`section9_residual_decay`.  This
module adds only an arithmetic witness layer: a caller must provide a genuinely
certified/theorem-derived upper bound and provenance for it, and this code
checks that bound against the corresponding paper envelope at one supplied
``q``.  It never samples a field, fits ``C/P/K``, manufactures ``E``, or turns a
pointwise check into a uniform-in-``q`` theorem.

The comparison is performed in logarithmic scale so very small powers of ``q``
do not underflow in binary64.  Passing a check therefore means only that the
supplied certified datum is arithmetically compatible with (9.17) or (9.18) at
that point.  Section 9 remains ``formal-structure`` until the genuine
correction sequence supplies these data uniformly on the common domain and the
flat remainder is proved to every order.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral

from .section9_residual_decay import (
    Scalar,
    section9_increment_derivative_loss,
    section9_increment_gain,
    section9_residual_gain,
)

_ALLOWED_EVIDENCE_KINDS = frozenset(
    {"paper-derived", "lean-derived", "certified-numerical", "rigorous-external"}
)


@dataclass(frozen=True)
class CertifiedBoundDatum:
    """One externally justified nonnegative upper bound.

    ``kind`` deliberately excludes ordinary sampled/fitted estimates.  The
    provenance string should identify the theorem, proof artifact, interval
    enclosure, or other source that justifies ``upper_bound``.
    """

    upper_bound: float
    kind: str
    provenance: str

    def __post_init__(self) -> None:
        upper = float(self.upper_bound)
        if not math.isfinite(upper) or upper < 0.0:
            raise ValueError("upper_bound must be finite and nonnegative")
        kind = str(self.kind).strip()
        if kind not in _ALLOWED_EVIDENCE_KINDS:
            raise ValueError(
                "kind must be one of " + ", ".join(sorted(_ALLOWED_EVIDENCE_KINDS))
            )
        provenance = str(self.provenance).strip()
        if not provenance:
            raise ValueError("provenance must be nonempty")
        object.__setattr__(self, "upper_bound", upper)
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "provenance", provenance)


@dataclass(frozen=True)
class Section9PointwiseBoundCertificate:
    """Arithmetic result for one supplied Eq. (9.17) or Eq. (9.18) datum."""

    equation: str
    stage: int
    derivative_order: int
    q: float
    exponent: Fraction
    constant: float
    log_power: int
    supplied_upper_bound: float
    leading_rhs_log: float
    total_rhs_log: float
    arithmetic_check_passed: bool
    evidence_kind: str
    evidence_provenance: str
    flat_remainder_upper_bound: float = 0.0
    flat_remainder_provenance: str | None = None
    status: str = "formal-structure"
    uniform_in_q_verified: bool = False
    flat_remainder_all_orders_verified: bool = False
    actual_section9_sequence_verified: bool = False
    paper_exact_velocity_available: bool = False


def _nonnegative_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _positive_constant(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _q(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or not 0.0 < value < 1.0:
        raise ValueError("q must be finite with 0 < q < 1")
    return value


def _as_fraction(value: Scalar, name: str) -> Fraction:
    """Convert theorem data exactly enough for audit arithmetic.

    Floats follow the decimal-string convention used by the existing Section 9
    exponent ledger rather than exposing their binary representation.  Boolean
    values are rejected explicitly instead of being accepted as integers.
    """

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real number")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return Fraction(str(value))


def _leading_log_rhs(
    q: float, exponent: Fraction, constant: float, log_power: int
) -> float:
    log_q = math.log(q)
    return (
        math.log(constant)
        + float(exponent) * log_q
        + log_power * math.log1p(abs(log_q))
    )


def _log_nonnegative(value: float) -> float:
    return -math.inf if value == 0.0 else math.log(value)


def _logaddexp(a: float, b: float) -> float:
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    hi = max(a, b)
    lo = min(a, b)
    return hi + math.log1p(math.exp(lo - hi))


def _compatible(upper_bound: float, rhs_log: float) -> bool:
    if upper_bound == 0.0:
        return True
    # The tolerance only protects the final binary64 comparison; it does not
    # enlarge the theorem constants or alter any exact exponent.
    tol = 32.0 * math.ulp(max(1.0, abs(rhs_log)))
    return math.log(upper_bound) <= rhs_log + tol


def check_eq_9_17_pointwise(
    *,
    stage: int,
    derivative_order: int,
    h: Scalar,
    q: float,
    constant: float,
    log_power: int,
    correction_bound: CertifiedBoundDatum,
) -> Section9PointwiseBoundCertificate:
    """Check one externally certified datum against Eq. (9.17) at ``q``.

    The supplied datum must bound ``|Z_j|_m + |Delta u_j|_m``.  This function
    computes no field norm itself and makes no uniform-in-``q`` claim.
    """

    stage = _nonnegative_integer(stage, "stage")
    derivative_order = _nonnegative_integer(derivative_order, "derivative_order")
    q = _q(q)
    constant = _positive_constant(constant, "constant")
    log_power = _nonnegative_integer(log_power, "log_power")
    if not isinstance(correction_bound, CertifiedBoundDatum):
        raise TypeError("correction_bound must be CertifiedBoundDatum")

    exponent = section9_increment_gain(stage, h) - section9_increment_derivative_loss(
        derivative_order, h
    )
    leading_log = _leading_log_rhs(q, exponent, constant, log_power)
    return Section9PointwiseBoundCertificate(
        equation="9.17",
        stage=stage,
        derivative_order=derivative_order,
        q=q,
        exponent=exponent,
        constant=constant,
        log_power=log_power,
        supplied_upper_bound=correction_bound.upper_bound,
        leading_rhs_log=leading_log,
        total_rhs_log=leading_log,
        arithmetic_check_passed=_compatible(correction_bound.upper_bound, leading_log),
        evidence_kind=correction_bound.kind,
        evidence_provenance=correction_bound.provenance,
    )


def check_eq_9_18_pointwise(
    *,
    stage: int,
    derivative_order: int,
    h: Scalar,
    q: float,
    derivative_loss: Scalar,
    constant: float,
    log_power: int,
    residual_bound: CertifiedBoundDatum,
    flat_remainder_bound: CertifiedBoundDatum,
) -> Section9PointwiseBoundCertificate:
    """Check one externally certified residual datum against Eq. (9.18) at ``q``.

    ``flat_remainder_bound`` is only a bound for the supplied point.  Even when
    this arithmetic check passes, the stronger paper statement
    ``E_{j,m} <= C_{j,m,N} q^N`` for every ``N`` remains unverified here.
    """

    stage = _nonnegative_integer(stage, "stage")
    derivative_order = _nonnegative_integer(derivative_order, "derivative_order")
    q = _q(q)
    constant = _positive_constant(constant, "constant")
    log_power = _nonnegative_integer(log_power, "log_power")
    if not isinstance(residual_bound, CertifiedBoundDatum):
        raise TypeError("residual_bound must be CertifiedBoundDatum")
    if not isinstance(flat_remainder_bound, CertifiedBoundDatum):
        raise TypeError("flat_remainder_bound must be CertifiedBoundDatum")

    loss_fraction = _as_fraction(derivative_loss, "derivative_loss")
    if loss_fraction < 0:
        raise ValueError("derivative_loss must be nonnegative")

    exponent = section9_residual_gain(stage, h) - loss_fraction
    leading_log = _leading_log_rhs(q, exponent, constant, log_power)
    flat_log = _log_nonnegative(flat_remainder_bound.upper_bound)
    total_log = _logaddexp(leading_log, flat_log)
    return Section9PointwiseBoundCertificate(
        equation="9.18",
        stage=stage,
        derivative_order=derivative_order,
        q=q,
        exponent=exponent,
        constant=constant,
        log_power=log_power,
        supplied_upper_bound=residual_bound.upper_bound,
        leading_rhs_log=leading_log,
        total_rhs_log=total_log,
        arithmetic_check_passed=_compatible(residual_bound.upper_bound, total_log),
        evidence_kind=residual_bound.kind,
        evidence_provenance=residual_bound.provenance,
        flat_remainder_upper_bound=flat_remainder_bound.upper_bound,
        flat_remainder_provenance=flat_remainder_bound.provenance,
    )
