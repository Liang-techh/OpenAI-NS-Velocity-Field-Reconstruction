"""Analytic composition bridge for the Section 9 cutoff weight ``chi(a_j q)``.

The weighted-correction adapter on main accepts a complete spacetime derivative
jet for ``chi(a_j q)``.  Supplying that final table directly leaves the chain
rule itself inside the provider trust boundary.  This module narrows that
boundary: callers provide

* a derivative jet for the scalar similarity coordinate ``q(t,x,y,z)``; and
* the one-variable derivatives ``chi^(k)(a_j q)`` at the cutoff argument.

The code then derives every spacetime derivative of ``chi(a_j q)`` through a
requested finite order by composing *normalized multivariate Taylor jets*.
This is an exact implementation of the multivariate chain/Faa-di-Bruno algebra;
no finite differences, fitting, or sampling are used.

This is deliberately only an analytic adapter.  The ``q`` derivatives and the
one-variable cutoff derivatives are still provider inputs.  In particular this
module does not prove that the q-jet came from Eq. (4.1), does not identify the
caller's scalar cutoff with the manuscript's fixed ``chi``, and does not upgrade
any endpoint/residual/paper-exact truth flag.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral, Rational, Real
from types import MappingProxyType
from typing import Mapping

from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_eq921_prefix_jet import MultiIndex4, required_spacetime_multiindices
from .section9_finite_prefix import Section9FinitePrefixEvaluationCertificate
from .section9_weighted_correction_jet import Section9CutoffWeightJet


_ZERO: MultiIndex4 = (0, 0, 0, 0)


def _fraction(value: object, name: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive finite rational value")
    if isinstance(value, Fraction):
        out = value
    elif isinstance(value, Integral):
        out = Fraction(int(value), 1)
    elif isinstance(value, Rational):
        out = Fraction(value)
    elif isinstance(value, Real):
        x = float(value)
        if not math.isfinite(x):
            raise ValueError(f"{name} must be finite")
        out = Fraction.from_float(x)
    else:
        raise TypeError(f"{name} must be a real rationalizable value")
    if out <= 0:
        raise ValueError(f"{name} must be positive")
    return out


def _nonnegative_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _nonempty(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty")
    return value.strip()


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _normalize_multiindex(raw: object, *, order: int, name: str) -> MultiIndex4:
    if not isinstance(raw, tuple) or len(raw) != 4:
        raise ValueError(f"{name} keys must be four-tuples (dt, dx, dy, dz)")
    if any(isinstance(v, bool) or not isinstance(v, Integral) for v in raw):
        raise ValueError(f"{name} keys must contain nonnegative integers")
    alpha = tuple(int(v) for v in raw)
    if any(v < 0 for v in alpha) or sum(alpha) > order:
        raise ValueError(f"{name} contains an invalid derivative multi-index")
    return alpha  # type: ignore[return-value]


def _normalize_q_derivatives(
    values: Mapping[MultiIndex4, object], *, order: int
) -> Mapping[MultiIndex4, float]:
    if not isinstance(values, Mapping):
        raise TypeError("q_derivatives must be a mapping")
    normalized: dict[MultiIndex4, float] = {}
    for raw_key, raw_value in values.items():
        key = _normalize_multiindex(raw_key, order=order, name="q_derivatives")
        if key in normalized:
            raise ValueError("q_derivatives contains a duplicate normalized multi-index")
        normalized[key] = _finite(raw_value, f"q_derivatives[{key}]")
    required = set(required_spacetime_multiindices(order))
    supplied = set(normalized)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            "q_derivatives must cover exactly every spacetime derivative through "
            f"total order {order}; missing={missing}, extra={extra}"
        )
    return MappingProxyType(dict(sorted(normalized.items())))


def _normalize_cutoff_derivatives(
    values: Mapping[int, object], *, order: int
) -> Mapping[int, float]:
    if not isinstance(values, Mapping):
        raise TypeError("cutoff_derivatives must be a mapping")
    normalized: dict[int, float] = {}
    for raw_key, raw_value in values.items():
        if isinstance(raw_key, bool) or not isinstance(raw_key, Integral):
            raise ValueError("cutoff_derivatives keys must be integer derivative orders")
        key = int(raw_key)
        if key < 0 or key > order:
            raise ValueError("cutoff_derivatives contains an order outside derivative_order")
        if key in normalized:
            raise ValueError("cutoff_derivatives contains a duplicate derivative order")
        normalized[key] = _finite(raw_value, f"cutoff_derivatives[{key}]")
    required = set(range(order + 1))
    supplied = set(normalized)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            "cutoff_derivatives must cover exactly orders 0..derivative_order; "
            f"missing={missing}, extra={extra}"
        )
    return MappingProxyType(dict(sorted(normalized.items())))


def _multi_factorial(alpha: MultiIndex4) -> int:
    return math.prod(math.factorial(component) for component in alpha)


def _add_multiindices(alpha: MultiIndex4, beta: MultiIndex4) -> MultiIndex4:
    return tuple(a + b for a, b in zip(alpha, beta))  # type: ignore[return-value]


def _multiply_normalized_jets(
    left: Mapping[MultiIndex4, float],
    right: Mapping[MultiIndex4, float],
    *,
    order: int,
) -> dict[MultiIndex4, float]:
    """Multiply normalized Taylor coefficients, truncating by total order."""
    out: dict[MultiIndex4, float] = {}
    for alpha, left_value in left.items():
        for beta, right_value in right.items():
            gamma = _add_multiindices(alpha, beta)
            if sum(gamma) <= order:
                out[gamma] = out.get(gamma, 0.0) + left_value * right_value
    if any(not math.isfinite(value) for value in out.values()):
        raise ArithmeticError("cutoff composition Taylor product overflowed")
    return out


@dataclass(frozen=True)
class Section9SimilarityCoordinateJet:
    """Provider input for a finite spacetime jet of ``q(t,x,y,z)``."""

    q: Fraction | float | int
    derivative_order: int
    derivatives: Mapping[MultiIndex4, object]
    provider_id: str
    provider_revision: str
    provider_provenance: str

    def __post_init__(self) -> None:
        q = _fraction(self.q, "q")
        order = _nonnegative_int(self.derivative_order, "derivative_order")
        derivatives = _normalize_q_derivatives(self.derivatives, order=order)
        if derivatives[_ZERO] != float(q):
            raise ValueError("q_derivatives zero-order value must equal q")
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "derivative_order", order)
        object.__setattr__(self, "derivatives", derivatives)
        for name in ("provider_id", "provider_revision", "provider_provenance"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), name))


@dataclass(frozen=True)
class Section9ScalarCutoffDerivativeJet:
    """Provider input for ``chi^(k)(a_j q)`` at one scalar argument."""

    argument: Fraction | float | int
    derivative_order: int
    derivatives: Mapping[int, object]
    provider_id: str
    provider_revision: str
    provider_provenance: str

    def __post_init__(self) -> None:
        argument = _fraction(self.argument, "argument")
        order = _nonnegative_int(self.derivative_order, "derivative_order")
        object.__setattr__(self, "argument", argument)
        object.__setattr__(self, "derivative_order", order)
        object.__setattr__(
            self,
            "derivatives",
            _normalize_cutoff_derivatives(self.derivatives, order=order),
        )
        for name in ("provider_id", "provider_revision", "provider_provenance"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), name))


@dataclass(frozen=True)
class Section9DerivedCutoffWeightJetCertificate:
    """Composition certificate plus the existing cutoff-weight jet payload."""

    weight_jet: Section9CutoffWeightJet
    q_provider_id: str
    q_provider_revision: str
    q_provider_provenance: str
    cutoff_provider_id: str
    cutoff_provider_revision: str
    cutoff_provider_provenance: str
    normalized_taylor_composition_verified: bool = True
    similarity_coordinate_jet_machine_derived_from_eq_4_1: bool = False
    paper_fixed_cutoff_derivatives_machine_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_composition_ready(self) -> bool:
        return (
            self.normalized_taylor_composition_verified
            and not self.similarity_coordinate_jet_machine_derived_from_eq_4_1
            and not self.paper_fixed_cutoff_derivatives_machine_verified
            and not self.paper_exact_velocity_available
        )


def derive_section9_cutoff_weight_jet(
    evaluation: Section9FinitePrefixEvaluationCertificate,
    admission: Section9CorrectionStageAdmissionCertificate,
    q_jet: Section9SimilarityCoordinateJet,
    cutoff_jet: Section9ScalarCutoffDerivativeJet,
) -> Section9DerivedCutoffWeightJetCertificate:
    """Derive all ``D^alpha chi(a_j q)`` through the supplied finite order.

    Normalized Taylor coefficients are ``D^alpha f / alpha!``.  Writing
    ``delta = a_j (q-q0)``, the finite jet is obtained from

        chi(a_j q0 + delta)
        = sum_k chi^(k)(a_j q0) delta^k / k!.

    Polynomial multiplication is truncated only above ``derivative_order``.
    This is algebraically equivalent to multivariate Faa di Bruno while being
    auditable and independent of finite-difference step sizes.
    """
    if not isinstance(evaluation, Section9FinitePrefixEvaluationCertificate):
        raise TypeError("evaluation must be a Section9FinitePrefixEvaluationCertificate")
    if not evaluation.formal_prefix_ready:
        raise ValueError("evaluation must be a complete formal finite-prefix evaluation")
    if not isinstance(admission, Section9CorrectionStageAdmissionCertificate):
        raise TypeError("admission must be a Section9CorrectionStageAdmissionCertificate")
    if not admission.formal_admission_ready:
        raise ValueError("admission must be a complete formal stage admission")
    if not isinstance(q_jet, Section9SimilarityCoordinateJet):
        raise TypeError("q_jet must be a Section9SimilarityCoordinateJet")
    if not isinstance(cutoff_jet, Section9ScalarCutoffDerivativeJet):
        raise TypeError("cutoff_jet must be a Section9ScalarCutoffDerivativeJet")

    stage = _positive_int(admission.stage, "stage")
    if q_jet.q != evaluation.q:
        raise ValueError("q jet must use the finite-prefix evaluation q")
    if q_jet.derivative_order != cutoff_jet.derivative_order:
        raise ValueError("q and scalar cutoff jets must use one derivative_order")

    rows = [row for row in evaluation.contributions if row.stage == stage]
    if len(rows) != 1:
        raise ValueError("evaluation must contain exactly one contribution for the admitted stage")
    contribution = rows[0]
    if contribution.scale != admission.scale:
        raise ValueError("evaluated contribution scale must match the admitted a_j")

    argument = admission.scale * evaluation.q
    if cutoff_jet.argument != argument:
        raise ValueError("scalar cutoff jet argument must equal the admitted a_j q")
    if cutoff_jet.provider_provenance != evaluation.cutoff_evaluator_provenance:
        raise ValueError("scalar cutoff provenance must match the finite-prefix cutoff evaluator")
    if cutoff_jet.derivatives[0] != contribution.cutoff_weight:
        raise ValueError("scalar cutoff zero-order value must match the finite-prefix cutoff weight")

    order = q_jet.derivative_order
    required = required_spacetime_multiindices(order)

    # Normalized jet for delta = a_j(q-q0).  The constant coefficient is zero.
    delta: dict[MultiIndex4, float] = {_ZERO: 0.0}
    for alpha in required:
        if alpha != _ZERO:
            delta[alpha] = (
                admission.scale * q_jet.derivatives[alpha] / _multi_factorial(alpha)
            )

    normalized_result = {alpha: 0.0 for alpha in required}
    power: dict[MultiIndex4, float] = {_ZERO: 1.0}
    for k in range(order + 1):
        coefficient = cutoff_jet.derivatives[k] / math.factorial(k)
        for alpha, value in power.items():
            normalized_result[alpha] += coefficient * value
        if k != order:
            power = _multiply_normalized_jets(power, delta, order=order)

    derivatives: dict[MultiIndex4, float] = {}
    for alpha in required:
        value = normalized_result[alpha] * _multi_factorial(alpha)
        if not math.isfinite(value):
            raise ArithmeticError("derived cutoff-weight derivative overflowed")
        derivatives[alpha] = float(value)

    weight_jet = Section9CutoffWeightJet(
        stage=stage,
        scale=admission.scale,
        q=evaluation.q,
        derivative_order=order,
        derivatives=derivatives,
        provider_id=(
            "analytic-composition:"
            f"{q_jet.provider_id}+{cutoff_jet.provider_id}"
        ),
        provider_revision=(
            f"q={q_jet.provider_revision};chi={cutoff_jet.provider_revision}"
        ),
        # Existing weighted-jet plumbing binds this field to the zero-order
        # cutoff evaluator, so preserve that evaluator's provenance verbatim.
        provider_provenance=cutoff_jet.provider_provenance,
    )
    certificate = Section9DerivedCutoffWeightJetCertificate(
        weight_jet=weight_jet,
        q_provider_id=q_jet.provider_id,
        q_provider_revision=q_jet.provider_revision,
        q_provider_provenance=q_jet.provider_provenance,
        cutoff_provider_id=cutoff_jet.provider_id,
        cutoff_provider_revision=cutoff_jet.provider_revision,
        cutoff_provider_provenance=cutoff_jet.provider_provenance,
    )
    if not certificate.formal_composition_ready:
        raise ArithmeticError("Section 9 cutoff composition truth-status invariant failed")
    return certificate
