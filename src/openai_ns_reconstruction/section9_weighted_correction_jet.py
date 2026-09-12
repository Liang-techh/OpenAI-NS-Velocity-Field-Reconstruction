"""Analytic product-rule bridge for one Section 9 Eq. (9.21) stage.

The finite-prefix jet adapter on main accepts derivatives of the already-weighted
summand ``chi(a_j q) (A_j, B_j, p_j)``.  That is useful for aggregation but
leaves a large verification gap: a provider could hand over the final weighted
derivative table directly.

This module narrows that gap by separating two inputs for one already-admitted
stage:

* an unweighted analytic spacetime jet for ``(A_j, B_j, p_j)``; and
* an analytic spacetime jet for the scalar cutoff weight ``chi(a_j q)``.

It then derives every derivative of the weighted summand with the exact
multi-index Leibniz rule and returns the existing
:class:`Section9WeightedCorrectionStageJet`.  The zero-order cutoff value and
weighted field values are cross-checked against the independently landed
finite-prefix evaluation before the result is accepted.

The input derivative tables remain provider data.  In particular, this module
does not derive the correction jets from the actual Section 7/8 construction,
does not derive the cutoff jet from the manuscript's implicit similarity
coordinate, and does not certify the manuscript's fixed cutoff.  Consequently
it does not promote any paper-exact truth flag or make a residual artifact.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import math
from numbers import Integral, Rational, Real
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_eq921_prefix_jet import (
    MultiIndex4,
    Section9WeightedCorrectionStageJet,
    required_spacetime_multiindices,
)
from .section9_finite_prefix import Section9FinitePrefixEvaluationCertificate


_ZERO_MULTIINDEX: MultiIndex4 = (0, 0, 0, 0)


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _nonnegative_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


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


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty")
    return value.strip()


def _multiindex(value: object, *, order: int, name: str) -> MultiIndex4:
    if not isinstance(value, tuple) or len(value) != 4:
        raise ValueError(f"{name} keys must be four-tuples (dt, dx, dy, dz)")
    if any(isinstance(item, bool) or not isinstance(item, Integral) for item in value):
        raise ValueError(f"{name} keys must contain nonnegative integers")
    alpha = tuple(int(item) for item in value)
    if any(item < 0 for item in alpha):
        raise ValueError(f"{name} keys must contain nonnegative integers")
    if sum(alpha) > order:
        raise ValueError(f"{name} contains a derivative above derivative_order")
    return alpha  # type: ignore[return-value]


def _finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    out = np.array(out, dtype=float, copy=True)
    out.setflags(write=False)
    return out


def _normalize_scalar_jet(
    values: Mapping[MultiIndex4, object], *, order: int, name: str
) -> Mapping[MultiIndex4, float]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{name} must be a mapping from spacetime multi-indices")
    normalized: dict[MultiIndex4, float] = {}
    for raw_key, raw_value in values.items():
        key = _multiindex(raw_key, order=order, name=name)
        if key in normalized:
            raise ValueError(f"{name} contains a duplicate normalized multi-index")
        normalized[key] = _finite_scalar(raw_value, f"{name}[{key}]")
    required = set(required_spacetime_multiindices(order))
    supplied = set(normalized)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            f"{name} must cover exactly all derivatives through total order {order}; "
            f"missing={missing}, extra={extra}"
        )
    return MappingProxyType(dict(sorted(normalized.items())))


def _normalize_vector_jet(
    values: Mapping[MultiIndex4, object], *, order: int, name: str
) -> Mapping[MultiIndex4, np.ndarray]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{name} must be a mapping from spacetime multi-indices")
    normalized: dict[MultiIndex4, np.ndarray] = {}
    for raw_key, raw_value in values.items():
        key = _multiindex(raw_key, order=order, name=name)
        if key in normalized:
            raise ValueError(f"{name} contains a duplicate normalized multi-index")
        normalized[key] = _vec3(raw_value, f"{name}[{key}]")
    required = set(required_spacetime_multiindices(order))
    supplied = set(normalized)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            f"{name} must cover exactly all derivatives through total order {order}; "
            f"missing={missing}, extra={extra}"
        )
    return MappingProxyType(dict(sorted(normalized.items())))


@dataclass(frozen=True)
class Section9UnweightedCorrectionStageJet:
    """Provider input for one unweighted ``(A_j,B_j,p_j)`` spacetime jet."""

    stage: int
    scale: int
    q: Fraction | float | int
    derivative_order: int
    A_derivatives: Mapping[MultiIndex4, object]
    B_derivatives: Mapping[MultiIndex4, object]
    p_derivatives: Mapping[MultiIndex4, object]
    evidence_kinds: tuple[str, str, str]
    evidence_provenance: tuple[str, str, str]
    value_provenance: str
    provider_id: str
    provider_revision: str
    provider_provenance: str

    def __post_init__(self) -> None:
        stage = _positive_int(self.stage, "stage")
        scale = _positive_int(self.scale, "scale")
        q = _fraction(self.q, "q")
        order = _nonnegative_int(self.derivative_order, "derivative_order")
        if (
            not isinstance(self.evidence_kinds, tuple)
            or len(self.evidence_kinds) != 3
            or any(not isinstance(value, str) or not value for value in self.evidence_kinds)
        ):
            raise ValueError("evidence_kinds must be a nonempty (A,B,p) string triple")
        if (
            not isinstance(self.evidence_provenance, tuple)
            or len(self.evidence_provenance) != 3
            or any(not isinstance(value, str) or not value for value in self.evidence_provenance)
        ):
            raise ValueError("evidence_provenance must be a nonempty (A,B,p) string triple")

        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "scale", scale)
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "derivative_order", order)
        object.__setattr__(
            self,
            "A_derivatives",
            _normalize_vector_jet(self.A_derivatives, order=order, name="A_derivatives"),
        )
        object.__setattr__(
            self,
            "B_derivatives",
            _normalize_scalar_jet(self.B_derivatives, order=order, name="B_derivatives"),
        )
        object.__setattr__(
            self,
            "p_derivatives",
            _normalize_scalar_jet(self.p_derivatives, order=order, name="p_derivatives"),
        )
        for name in (
            "value_provenance",
            "provider_id",
            "provider_revision",
            "provider_provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))


@dataclass(frozen=True)
class Section9CutoffWeightJet:
    """Provider input for derivatives of the scalar weight ``chi(a_j q)``."""

    stage: int
    scale: int
    q: Fraction | float | int
    derivative_order: int
    derivatives: Mapping[MultiIndex4, object]
    provider_id: str
    provider_revision: str
    provider_provenance: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "stage", _positive_int(self.stage, "stage"))
        object.__setattr__(self, "scale", _positive_int(self.scale, "scale"))
        object.__setattr__(self, "q", _fraction(self.q, "q"))
        order = _nonnegative_int(self.derivative_order, "derivative_order")
        object.__setattr__(self, "derivative_order", order)
        object.__setattr__(
            self,
            "derivatives",
            _normalize_scalar_jet(self.derivatives, order=order, name="cutoff_derivatives"),
        )
        for name in ("provider_id", "provider_revision", "provider_provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))


def _sub_multiindices(alpha: MultiIndex4) -> tuple[MultiIndex4, ...]:
    return tuple(
        tuple(int(value) for value in beta)  # type: ignore[misc]
        for beta in product(*(range(component + 1) for component in alpha))
    )


def _multi_binomial(alpha: MultiIndex4, beta: MultiIndex4) -> int:
    return math.prod(math.comb(a, b) for a, b in zip(alpha, beta))


def _difference(alpha: MultiIndex4, beta: MultiIndex4) -> MultiIndex4:
    return tuple(a - b for a, b in zip(alpha, beta))  # type: ignore[return-value]


def _weighted_scalar_derivative(
    alpha: MultiIndex4,
    weight: Mapping[MultiIndex4, float],
    field: Mapping[MultiIndex4, float],
) -> float:
    out = 0.0
    for beta in _sub_multiindices(alpha):
        gamma = _difference(alpha, beta)
        out += _multi_binomial(alpha, beta) * weight[beta] * field[gamma]
    if not math.isfinite(out):
        raise ArithmeticError("weighted scalar derivative overflowed")
    return float(out)


def _weighted_vector_derivative(
    alpha: MultiIndex4,
    weight: Mapping[MultiIndex4, float],
    field: Mapping[MultiIndex4, np.ndarray],
) -> np.ndarray:
    out = np.zeros(3, dtype=float)
    for beta in _sub_multiindices(alpha):
        gamma = _difference(alpha, beta)
        out += _multi_binomial(alpha, beta) * weight[beta] * field[gamma]
    if not np.all(np.isfinite(out)):
        raise ArithmeticError("weighted vector derivative overflowed")
    out.setflags(write=False)
    return out


def derive_weighted_correction_stage_jet(
    evaluation: Section9FinitePrefixEvaluationCertificate,
    admission: Section9CorrectionStageAdmissionCertificate,
    correction_jet: Section9UnweightedCorrectionStageJet,
    cutoff_jet: Section9CutoffWeightJet,
) -> Section9WeightedCorrectionStageJet:
    """Derive ``D^alpha[chi(a_j q) (A_j,B_j,p_j)]`` by Leibniz rule.

    The returned weighted jet can be passed directly to
    :func:`assemble_eq_9_21_prefix_jet`.  This function verifies only analytic
    product-rule consistency and source identity.  It deliberately does not
    certify that either input jet came from the paper's actual construction.
    """
    if not isinstance(evaluation, Section9FinitePrefixEvaluationCertificate):
        raise TypeError("evaluation must be a Section9FinitePrefixEvaluationCertificate")
    if not evaluation.formal_prefix_ready:
        raise ValueError("evaluation must be a complete formal finite-prefix evaluation")
    if not isinstance(admission, Section9CorrectionStageAdmissionCertificate):
        raise TypeError("admission must be a Section9CorrectionStageAdmissionCertificate")
    if not admission.formal_admission_ready:
        raise ValueError("admission must be a complete formal stage admission")
    if not isinstance(correction_jet, Section9UnweightedCorrectionStageJet):
        raise TypeError("correction_jet must be a Section9UnweightedCorrectionStageJet")
    if not isinstance(cutoff_jet, Section9CutoffWeightJet):
        raise TypeError("cutoff_jet must be a Section9CutoffWeightJet")

    stage = admission.stage
    if correction_jet.stage != stage or cutoff_jet.stage != stage:
        raise ValueError("admission, correction jet and cutoff jet must use one stage")
    if correction_jet.scale != admission.scale or cutoff_jet.scale != admission.scale:
        raise ValueError("correction and cutoff jet scales must match the admitted a_j")
    if correction_jet.q != evaluation.q or cutoff_jet.q != evaluation.q:
        raise ValueError("correction and cutoff jets must use the finite-prefix evaluation q")
    if correction_jet.derivative_order != cutoff_jet.derivative_order:
        raise ValueError("correction and cutoff jets must use one derivative_order")
    if correction_jet.evidence_kinds != admission.evidence_kinds:
        raise ValueError("correction jet evidence kinds must match the stage admission")
    if correction_jet.evidence_provenance != admission.evidence_provenance:
        raise ValueError("correction jet evidence provenance must match the stage admission")

    contributions = [row for row in evaluation.contributions if row.stage == stage]
    if len(contributions) != 1:
        raise ValueError("evaluation must contain exactly one contribution for the admitted stage")
    contribution = contributions[0]
    if contribution.scale != admission.scale:
        raise ValueError("evaluated contribution scale must match the admitted a_j")
    if contribution.value_provenance != correction_jet.value_provenance:
        raise ValueError("correction jet value provenance must match the finite-prefix contribution")
    if evaluation.cutoff_evaluator_provenance != cutoff_jet.provider_provenance:
        raise ValueError("cutoff jet provenance must match the finite-prefix cutoff evaluator")
    if cutoff_jet.derivatives[_ZERO_MULTIINDEX] != contribution.cutoff_weight:
        raise ValueError("cutoff jet zero-order value must match the finite-prefix cutoff weight")

    order = correction_jet.derivative_order
    A_weighted: dict[MultiIndex4, np.ndarray] = {}
    B_weighted: dict[MultiIndex4, float] = {}
    p_weighted: dict[MultiIndex4, float] = {}
    for alpha in required_spacetime_multiindices(order):
        A_weighted[alpha] = _weighted_vector_derivative(
            alpha, cutoff_jet.derivatives, correction_jet.A_derivatives
        )
        B_weighted[alpha] = _weighted_scalar_derivative(
            alpha, cutoff_jet.derivatives, correction_jet.B_derivatives
        )
        p_weighted[alpha] = _weighted_scalar_derivative(
            alpha, cutoff_jet.derivatives, correction_jet.p_derivatives
        )

    if not np.array_equal(A_weighted[_ZERO_MULTIINDEX], contribution.A_term):
        raise ValueError("product-rule zero-order A term does not match finite-prefix evaluation")
    if B_weighted[_ZERO_MULTIINDEX] != contribution.B_term:
        raise ValueError("product-rule zero-order B term does not match finite-prefix evaluation")
    if p_weighted[_ZERO_MULTIINDEX] != contribution.p_term:
        raise ValueError("product-rule zero-order p term does not match finite-prefix evaluation")

    provider_id = (
        "analytic-leibniz["
        f"correction={correction_jet.provider_id};cutoff={cutoff_jet.provider_id}]"
    )
    provider_revision = (
        f"correction={correction_jet.provider_revision};cutoff={cutoff_jet.provider_revision}"
    )
    provider_provenance = (
        "machine-derived multi-index Leibniz product from unweighted correction jet "
        f"({correction_jet.provider_provenance}) and cutoff-weight jet "
        f"({cutoff_jet.provider_provenance}); inputs remain provider-supplied"
    )

    return Section9WeightedCorrectionStageJet(
        stage=stage,
        scale=admission.scale,
        q=evaluation.q,
        derivative_order=order,
        A_derivatives=A_weighted,
        B_derivatives=B_weighted,
        p_derivatives=p_weighted,
        evidence_kinds=admission.evidence_kinds,
        evidence_provenance=admission.evidence_provenance,
        provider_id=provider_id,
        provider_revision=provider_revision,
        provider_provenance=provider_provenance,
    )
