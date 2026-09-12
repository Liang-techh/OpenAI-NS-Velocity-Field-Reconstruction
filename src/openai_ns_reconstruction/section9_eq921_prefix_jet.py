"""Analytic-jet plumbing for one admitted finite prefix of Section 9 Eq. (9.21).

The existing finite-prefix evaluator carries only point values of the weighted
summands ``chi(a_j q) (A_j, B_j, p_j)``. That is insufficient for a genuine
Navier--Stokes residual, which needs spacetime derivatives of the constructed
fields. This module adds a deliberately narrow adapter for such derivatives.

A caller may supply a *complete analytic derivative jet* for every weighted
summand in an already-admitted contiguous prefix. The adapter checks exact
stage/scale/q identity, exact admission provenance, one provider revision,
complete four-variable multi-index coverage, and exact agreement of zero-order
weighted summands with the independently landed finite-prefix evaluation. It
then adds the derivative jets term-by-term. No finite differences, fits, or
sampled derivative reconstructions are used here.

This remains an adapter, not the missing paper construction: the supplied jets
are not machine-derived here from the actual Section 7/8 correction fields or
from the manuscript cutoff. Therefore no infinite Eq. (9.21) sum, residual
artifact, endpoint smoothness, or paper-exact velocity is claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Rational, Real
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_finite_prefix import Section9FinitePrefixEvaluationCertificate


MultiIndex4 = tuple[int, int, int, int]
_ZERO_MULTIINDEX: MultiIndex4 = (0, 0, 0, 0)
_WEIGHTED_SUMMAND_DEFINITION = "chi(a_j q)*(A_j,B_j,p_j)"


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


def required_spacetime_multiindices(order: int) -> tuple[MultiIndex4, ...]:
    """Return every ``(dt, dx, dy, dz)`` with total degree at most ``order``."""
    degree = _nonnegative_int(order, "derivative_order")
    out: list[MultiIndex4] = []
    for dt in range(degree + 1):
        for dx in range(degree - dt + 1):
            for dy in range(degree - dt - dx + 1):
                for dz in range(degree - dt - dx - dy + 1):
                    out.append((dt, dx, dy, dz))
    return tuple(sorted(out, key=lambda alpha: (sum(alpha), alpha)))


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


def _vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    out = np.array(out, dtype=float, copy=True)
    out.setflags(write=False)
    return out


def _finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


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


@dataclass(frozen=True)
class Section9WeightedCorrectionStageJet:
    """Provider-supplied jet of one weighted Eq. (9.21) positive-stage summand.

    Derivatives are physical ``(t,x,y,z)`` derivatives. The maps must contain
    every multi-index through ``derivative_order``. Zero-order entries are later
    cross-checked against the existing finite-prefix evaluator; higher
    derivatives remain provider data.
    """

    stage: int
    scale: int
    q: Fraction | float | int
    derivative_order: int
    A_derivatives: Mapping[MultiIndex4, object]
    B_derivatives: Mapping[MultiIndex4, object]
    p_derivatives: Mapping[MultiIndex4, object]
    evidence_kinds: tuple[str, str, str]
    evidence_provenance: tuple[str, str, str]
    provider_id: str
    provider_revision: str
    provider_provenance: str
    summand_definition: str = _WEIGHTED_SUMMAND_DEFINITION

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
        if self.summand_definition != _WEIGHTED_SUMMAND_DEFINITION:
            raise ValueError("summand_definition must match the pinned Eq. (9.21) weighted term")

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
        for name in ("provider_id", "provider_revision", "provider_provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))


@dataclass(frozen=True)
class Section9FinitePrefixJetCertificate:
    """Summed derivative jet for one admitted finite Eq. (9.21) prefix."""

    q: Fraction
    prefix_order: int
    stages: tuple[int, ...]
    derivative_order: int
    A_derivatives: Mapping[MultiIndex4, np.ndarray]
    B_derivatives: Mapping[MultiIndex4, float]
    p_derivatives: Mapping[MultiIndex4, float]
    provider_id: str
    provider_revision: str
    provider_provenance: str
    complete_multiindex_coverage_verified: bool = True
    stage_scale_identity_verified: bool = True
    stage_evidence_identity_verified: bool = True
    single_provider_revision_verified: bool = True
    zero_order_matches_finite_prefix_evaluation: bool = True
    analytic_derivatives_machine_derived_from_actual_corrections: bool = False
    paper_fixed_cutoff_derivatives_machine_verified: bool = False
    actual_correction_field_values_verified: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    endpoint_covered: bool = False
    residual_artifact_ready: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_prefix_jet_ready(self) -> bool:
        required = set(required_spacetime_multiindices(self.derivative_order))
        return (
            self.prefix_order == len(self.stages)
            and self.stages == tuple(range(1, self.prefix_order + 1))
            and set(self.A_derivatives) == required
            and set(self.B_derivatives) == required
            and set(self.p_derivatives) == required
            and self.complete_multiindex_coverage_verified
            and self.stage_scale_identity_verified
            and self.stage_evidence_identity_verified
            and self.single_provider_revision_verified
            and self.zero_order_matches_finite_prefix_evaluation
            and not self.analytic_derivatives_machine_derived_from_actual_corrections
            and not self.paper_fixed_cutoff_derivatives_machine_verified
            and not self.actual_correction_field_values_verified
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.endpoint_covered
            and not self.residual_artifact_ready
            and not self.paper_exact_velocity_available
        )


def assemble_eq_9_21_prefix_jet(
    evaluation: Section9FinitePrefixEvaluationCertificate,
    stage_certificates: tuple[Section9CorrectionStageAdmissionCertificate, ...],
    stage_jets: tuple[Section9WeightedCorrectionStageJet, ...],
) -> Section9FinitePrefixJetCertificate:
    """Add analytic jets for exactly the weighted summands in ``evaluation``.

    Each supplied zero-order stage term must equal the corresponding independent
    finite-prefix contribution. Higher derivatives are not independently proved
    here; they remain fail-closed provider inputs until a genuine Section 7/8
    correction-field exporter supplies them.
    """
    if not isinstance(evaluation, Section9FinitePrefixEvaluationCertificate):
        raise TypeError("evaluation must be a Section9FinitePrefixEvaluationCertificate")
    if not evaluation.formal_prefix_ready:
        raise ValueError("evaluation must be a complete formal finite-prefix evaluation")

    cert_rows = tuple(stage_certificates)
    jet_rows = tuple(stage_jets)
    if not cert_rows or not jet_rows:
        raise ValueError("stage certificates and stage jets must be nonempty")
    if not all(isinstance(row, Section9CorrectionStageAdmissionCertificate) for row in cert_rows):
        raise TypeError("every stage certificate must be a Section9CorrectionStageAdmissionCertificate")
    if not all(isinstance(row, Section9WeightedCorrectionStageJet) for row in jet_rows):
        raise TypeError("every stage jet must be a Section9WeightedCorrectionStageJet")

    cert_by_stage: dict[int, Section9CorrectionStageAdmissionCertificate] = {}
    for row in cert_rows:
        if row.stage in cert_by_stage:
            raise ValueError(f"duplicate admitted Section 9 stage {row.stage}")
        if not row.formal_admission_ready:
            raise ValueError(f"stage {row.stage} is not a complete formal admission")
        cert_by_stage[row.stage] = row

    jet_by_stage: dict[int, Section9WeightedCorrectionStageJet] = {}
    for row in jet_rows:
        if row.stage in jet_by_stage:
            raise ValueError(f"duplicate Section 9 stage jet {row.stage}")
        jet_by_stage[row.stage] = row

    expected = evaluation.stages
    if tuple(sorted(cert_by_stage)) != expected:
        raise ValueError("stage certificates must cover exactly the evaluated contiguous prefix")
    if tuple(sorted(jet_by_stage)) != expected:
        raise ValueError("stage jets must cover exactly the evaluated contiguous prefix")
    if tuple(row.stage for row in evaluation.contributions) != expected:
        raise ValueError("finite-prefix contributions must be ordered by the evaluated stages")

    first_jet = jet_by_stage[expected[0]]
    common_order = first_jet.derivative_order
    common_provider = (
        first_jet.provider_id,
        first_jet.provider_revision,
        first_jet.provider_provenance,
    )

    contributions = {row.stage: row for row in evaluation.contributions}
    for stage in expected:
        cert = cert_by_stage[stage]
        jet = jet_by_stage[stage]
        contribution = contributions[stage]
        if jet.q != evaluation.q:
            raise ValueError("every stage jet must use the finite-prefix evaluation q")
        if jet.scale != cert.scale or jet.scale != contribution.scale:
            raise ValueError("stage jet scale must match both admission and evaluated contribution")
        if jet.derivative_order != common_order:
            raise ValueError("all stage jets must use one common derivative_order")
        if (
            jet.provider_id,
            jet.provider_revision,
            jet.provider_provenance,
        ) != common_provider:
            raise ValueError("all stage jets must come from one provider id/revision/provenance")
        if jet.evidence_kinds != cert.evidence_kinds:
            raise ValueError("stage jet evidence kinds must match the stage admission exactly")
        if jet.evidence_provenance != cert.evidence_provenance:
            raise ValueError("stage jet evidence provenance must match the stage admission exactly")
        if not np.array_equal(jet.A_derivatives[_ZERO_MULTIINDEX], contribution.A_term):
            raise ValueError("stage jet zero-order A term does not match finite-prefix evaluation")
        if jet.B_derivatives[_ZERO_MULTIINDEX] != contribution.B_term:
            raise ValueError("stage jet zero-order B term does not match finite-prefix evaluation")
        if jet.p_derivatives[_ZERO_MULTIINDEX] != contribution.p_term:
            raise ValueError("stage jet zero-order p term does not match finite-prefix evaluation")

    A_sum: dict[MultiIndex4, np.ndarray] = {}
    B_sum: dict[MultiIndex4, float] = {}
    p_sum: dict[MultiIndex4, float] = {}
    for alpha in required_spacetime_multiindices(common_order):
        A_value = np.zeros(3, dtype=float)
        B_value = 0.0
        p_value = 0.0
        for stage in expected:
            jet = jet_by_stage[stage]
            A_value += jet.A_derivatives[alpha]
            B_value += jet.B_derivatives[alpha]
            p_value += jet.p_derivatives[alpha]
        if not (np.all(np.isfinite(A_value)) and math.isfinite(B_value) and math.isfinite(p_value)):
            raise ArithmeticError("finite-prefix derivative-jet accumulation overflowed")
        A_value.setflags(write=False)
        A_sum[alpha] = A_value
        B_sum[alpha] = float(B_value)
        p_sum[alpha] = float(p_value)

    result = Section9FinitePrefixJetCertificate(
        q=evaluation.q,
        prefix_order=evaluation.prefix_order,
        stages=evaluation.stages,
        derivative_order=common_order,
        A_derivatives=MappingProxyType(A_sum),
        B_derivatives=MappingProxyType(B_sum),
        p_derivatives=MappingProxyType(p_sum),
        provider_id=common_provider[0],
        provider_revision=common_provider[1],
        provider_provenance=common_provider[2],
    )
    if not result.formal_prefix_jet_ready:
        raise ArithmeticError("Section 9 finite-prefix jet truth-status invariant failed")
    return result
