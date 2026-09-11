"""Formal-structure scale schedule for the Section 10 spatial Borel extension.

The pinned ``SpatialBorelExtension.lean`` chooses, for Taylor degree ``j``, a
natural number strictly larger than

    2**j * sum_{m<j} sum_{k<j} templateBound(m, j, k)

and then applies ``DiagonalScale.doublingEnvelope``.  The resulting scale is
large enough to make every localized derivative with ``m<j`` and ``k<j`` at
most ``2**(-j)``.

This module makes that *scale-selection arithmetic* executable.  It does not
construct or certify the analytic ``templateBound`` values themselves.  A
caller-supplied bound provider is an explicit dependency boundary: until those
bounds are derived from the actual residual endpoint jets, this remains
formal-structure infrastructure and is not evidence of a paper-exact force or
of smooth gluing through ``t=1``.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
import math
import operator

from .endpoint_borel import BorelRightExtension, DoublingEnvelope
from .endpoint_jets import FullSpacetimeJetFamily, normal_jet_family

TemplateBoundProvider = Callable[[int, int, int], object]


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(out)


def _nonnegative_fraction(value: object, name: str) -> Fraction:
    """Validate a nonnegative real-like bound and retain exact binary64 value.

    Integer/Fraction inputs stay exact.  Other real-like inputs are first
    converted to a finite Python float and then represented by its exact binary
    rational.  This prevents the schedule inequalities themselves from being
    weakened by a second round of floating-point arithmetic.
    """
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite nonnegative bound")
    if isinstance(value, Fraction):
        out = value
    elif isinstance(value, int):
        out = Fraction(value)
    else:
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{name} must be a finite nonnegative bound") from exc
        if not math.isfinite(numeric):
            raise ValueError(f"{name} must be a finite nonnegative bound")
        out = Fraction.from_float(numeric)
    if out < 0:
        raise ValueError(f"{name} must be a finite nonnegative bound")
    return out


@dataclass(frozen=True)
class DegreeScaleCertificate:
    """Exact arithmetic certificate for one Taylor degree's scale choice."""

    degree: int
    bound_sum: Fraction
    local_scale: int
    scale: int

    @property
    def local_target(self) -> Fraction:
        return (1 << self.degree) * self.bound_sum

    @property
    def geometric_scale_budget(self) -> Fraction:
        # Rearranged form of boundSum <= 2^{-j} * scale.
        return Fraction(self.scale, 1 << self.degree)

    @property
    def certified(self) -> bool:
        return (
            self.degree >= 0
            and self.bound_sum >= 0
            and self.local_scale > self.local_target
            and self.scale >= self.local_scale
            and self.scale >= 1
            and self.bound_sum <= self.geometric_scale_budget
        )


@dataclass(frozen=True)
class TailScaleCertificate:
    """Exact supplied-bound version of ``localized_derivative_tail_bound``."""

    spatial_window: int
    degree: int
    derivative_order: int
    template_bound: Fraction
    scale: int
    direct_upper_bound: Fraction
    geometric_target: Fraction

    @property
    def certified(self) -> bool:
        return (
            0 <= self.spatial_window < self.degree
            and 0 <= self.derivative_order < self.degree
            and self.template_bound >= 0
            and self.scale >= 1
            and self.direct_upper_bound <= self.geometric_target
        )


class SpatialBorelScaleSchedule:
    """Deterministic constructive witness for the pinned Borel scale recurrence.

    ``template_bound(m,j,k)`` is intended to provide a valid global bound for
    the order-``k`` derivative of ``SpatialBorelExtension.template m j (a j)``.
    The Lean source obtains such numbers nonconstructively from compact support
    and smoothness.  This class checks and consumes supplied values; it does not
    infer them from point samples.

    Lean uses ``Classical.choose`` for any natural strictly above its target.
    Here ``floor(target)+1`` is a deterministic admissible witness, not a claim
    to reproduce the opaque chosen natural definitionally.
    """

    def __init__(self, template_bound: TemplateBoundProvider):
        if not callable(template_bound):
            raise ValueError("template_bound must be callable")
        self._provider = template_bound
        self._bound_cache: dict[tuple[int, int, int], Fraction] = {}
        self._sum_cache: dict[int, Fraction] = {}
        self._envelope = DoublingEnvelope(self.local_scale)

    def template_bound(self, m: int, j: int, k: int) -> Fraction:
        m = _natural(m, "spatial window")
        j = _natural(j, "degree")
        k = _natural(k, "derivative order")
        key = (m, j, k)
        if key not in self._bound_cache:
            self._bound_cache[key] = _nonnegative_fraction(
                self._provider(m, j, k), f"template_bound({m},{j},{k})"
            )
        return self._bound_cache[key]

    def bound_sum(self, degree: int) -> Fraction:
        """Exact ``sum_{m<j} sum_{k<j} templateBound(m,j,k)``."""
        degree = _natural(degree, "degree")
        if degree not in self._sum_cache:
            total = Fraction(0)
            for m in range(degree):
                for k in range(degree):
                    total += self.template_bound(m, degree, k)
            self._sum_cache[degree] = total
        return self._sum_cache[degree]

    def local_scale(self, degree: int) -> int:
        """Least natural strictly above ``2**j * boundSum(j)``."""
        degree = _natural(degree, "degree")
        target = (1 << degree) * self.bound_sum(degree)
        return target.numerator // target.denominator + 1

    def scale(self, degree: int) -> int:
        """Pinned ``DiagonalScale.doublingEnvelope(localScale)`` recurrence."""
        return self._envelope(_natural(degree, "degree"))

    def degree_certificate(self, degree: int) -> DegreeScaleCertificate:
        degree = _natural(degree, "degree")
        certificate = DegreeScaleCertificate(
            degree=degree,
            bound_sum=self.bound_sum(degree),
            local_scale=self.local_scale(degree),
            scale=self.scale(degree),
        )
        if not certificate.certified:
            raise ArithmeticError("Spatial Borel degree-scale invariant failed")
        return certificate

    def tail_certificate(self, m: int, degree: int, k: int) -> TailScaleCertificate:
        """Certify the geometric tail inequality for one ``m<j, k<j`` triple."""
        m = _natural(m, "spatial window")
        degree = _natural(degree, "degree")
        k = _natural(k, "derivative order")
        if m >= degree or k >= degree:
            raise ValueError("tail certificate requires spatial window < degree and derivative order < degree")
        bound = self.template_bound(m, degree, k)
        scale = self.scale(degree)
        # Since k<j, (scale^k / scale^j) * C = C / scale^(j-k).
        direct = bound / (scale ** (degree - k))
        target = Fraction(1, 1 << degree)
        certificate = TailScaleCertificate(
            spatial_window=m,
            degree=degree,
            derivative_order=k,
            template_bound=bound,
            scale=scale,
            direct_upper_bound=direct,
            geometric_target=target,
        )
        if not certificate.certified:
            raise ArithmeticError("supplied template bounds do not satisfy the pinned tail inequality")
        return certificate


def borel_extension_from_template_bounds(
    full_jets: FullSpacetimeJetFamily,
    template_bound: TemplateBoundProvider,
    *,
    endpoint: float = 1.0,
) -> tuple[BorelRightExtension, SpatialBorelScaleSchedule]:
    """Wire the pinned scale arithmetic into the existing endpoint evaluator.

    This helper still accepts caller-supplied full endpoint jets and analytic
    template bounds.  It therefore does not assert ``CandidateFromLimits``'s
    locally-uniform residual derivative limits or the smoothness hypotheses
    from which Lean obtains ``templateBound``.
    """
    schedule = SpatialBorelScaleSchedule(template_bound)
    extension = BorelRightExtension(
        normal_jet_family(full_jets), schedule.local_scale, endpoint=endpoint
    )
    return extension, schedule
