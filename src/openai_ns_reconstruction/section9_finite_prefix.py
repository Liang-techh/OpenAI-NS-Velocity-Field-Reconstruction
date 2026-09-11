"""Fail-closed finite-prefix evaluation for the positive-stage part of Eq. (9.21).

Proposition 9.9 ultimately uses locally finite *infinite* sums of the admitted
Section 9 corrections.  The already-landed admission modules certify structural
facts about one stage (common domain, actual schedule membership, cutoff support
and smooth zero extension), but deliberately do not inspect correction values or
construct the Eq. (9.21) fields.

This module adds only the next plumbing layer: at one point with ``0 < q < q_big``
it can evaluate a *finite contiguous prefix* of already-admitted positive stages,

    sum_{j=1}^N chi(a_j q) A_j,
    sum_{j=1}^N chi(a_j q) B_j,
    sum_{j=1}^N chi(a_j q) p_j.

The numerical values and cutoff evaluator remain external inputs.  Consequently
the result is formal-structure only: it does not verify that supplied values are
the paper's actual correction fields, does not certify the pointwise evaluator
as the paper's fixed ``chi``, and never promotes a finite prefix to the infinite
locally finite sum or Proposition 9.9.

For ``a_j q >= 1`` the paper cutoff-support theorem already carried by the stage
admission is enough to short-circuit that entire stage to exact zero without
calling the external cutoff evaluator.  For ``a_j q < 1`` the evaluator is used
only as a numerical value and must return a finite number in ``[0,1]``.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Rational, Real
import math
from typing import Callable, Iterable

import numpy as np

from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_local_sum_admission import Section9LocalSumHypothesis


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


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


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


@dataclass(frozen=True)
class Section9CorrectionStageValue:
    """Numerical value payload for one positive Section 9 correction stage.

    ``provenance`` records where the numerical payload came from, but this class
    does not prove that it is the field covered by the corresponding admission
    theorem.  That separation is intentional and is reflected by the returned
    certificate's fail-closed truth flags.
    """

    stage: int
    A: np.ndarray
    B: float
    p: float
    provenance: str

    def __post_init__(self) -> None:
        stage = _positive_int(self.stage, "stage")
        A = _vec3(self.A, "A")
        B = _finite_scalar(self.B, "B")
        p = _finite_scalar(self.p, "p")
        provenance = str(self.provenance).strip()
        if not provenance:
            raise ValueError("value provenance must be nonempty")
        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "provenance", provenance)


@dataclass(frozen=True)
class Section9FinitePrefixStageContribution:
    stage: int
    scale: int
    cutoff_argument: Fraction
    cutoff_weight: float
    A_term: np.ndarray
    B_term: float
    p_term: float
    value_provenance: str
    support_zero_short_circuit: bool


@dataclass(frozen=True)
class Section9FinitePrefixEvaluationCertificate:
    """A finite-prefix numerical evaluation, never an infinite-sum certificate."""

    q: Fraction
    q_big: Fraction
    prefix_order: int
    stages: tuple[int, ...]
    contributions: tuple[Section9FinitePrefixStageContribution, ...]
    A_correction_prefix: np.ndarray
    B_correction_prefix: float
    p_correction_prefix: float
    cutoff_evaluator_provenance: str
    contiguous_prefix_verified: bool = True
    admitted_stage_structure_verified: bool = True
    paper_cutoff_support_short_circuits_verified: bool = True
    finite_prefix_evaluated: bool = True
    actual_correction_field_values_verified: bool = False
    paper_fixed_cutoff_pointwise_evaluator_verified: bool = False
    infinite_schedule_constructed_here: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    proposition_9_9_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_prefix_ready(self) -> bool:
        return (
            self.prefix_order == len(self.stages) == len(self.contributions)
            and self.stages == tuple(range(1, self.prefix_order + 1))
            and self.contiguous_prefix_verified
            and self.admitted_stage_structure_verified
            and self.paper_cutoff_support_short_circuits_verified
            and self.finite_prefix_evaluated
            and not self.actual_correction_field_values_verified
            and not self.paper_fixed_cutoff_pointwise_evaluator_verified
            and not self.infinite_schedule_constructed_here
            and not self.eq_9_21_infinite_sum_constructed
            and not self.proposition_9_9_verified
            and not self.paper_exact_velocity_available
        )


def evaluate_admitted_eq_9_21_prefix(
    hypothesis: Section9LocalSumHypothesis,
    stage_certificates: Iterable[Section9CorrectionStageAdmissionCertificate],
    stage_values: Iterable[Section9CorrectionStageValue],
    *,
    q: Fraction | float | int,
    cutoff_evaluator: Callable[[float], float],
    cutoff_evaluator_provenance: str,
) -> Section9FinitePrefixEvaluationCertificate:
    """Evaluate a contiguous admitted positive-stage prefix at one ``q``.

    The function verifies structural coherence against ``hypothesis`` and pairs
    exactly one numerical value payload with every admitted stage ``1..N``.  It
    refuses sparse stage sets because a sparse finite sum must not be presented
    as a truncation prefix.

    ``cutoff_evaluator`` is deliberately an external numerical evaluator.  This
    module cannot infer from a Python callable that it is the manuscript's fixed
    smooth cutoff, so the returned truth flag for that fact remains false.
    """
    if not isinstance(hypothesis, Section9LocalSumHypothesis):
        raise TypeError("hypothesis must be a Section9LocalSumHypothesis")
    q_value = _fraction(q, "q")
    if not q_value < hypothesis.q_big:
        raise ValueError("finite-prefix evaluation requires 0 < q < q_big")
    if not callable(cutoff_evaluator):
        raise TypeError("cutoff_evaluator must be callable")
    cutoff_provenance = str(cutoff_evaluator_provenance).strip()
    if not cutoff_provenance:
        raise ValueError("cutoff_evaluator_provenance must be nonempty")

    cert_rows = tuple(stage_certificates)
    value_rows = tuple(stage_values)
    if not cert_rows:
        raise ValueError("finite Section 9 prefix must contain at least stage 1")
    if not all(isinstance(row, Section9CorrectionStageAdmissionCertificate) for row in cert_rows):
        raise TypeError("every stage certificate must be a Section9CorrectionStageAdmissionCertificate")
    if not all(isinstance(row, Section9CorrectionStageValue) for row in value_rows):
        raise TypeError("every stage value must be a Section9CorrectionStageValue")

    cert_by_stage: dict[int, Section9CorrectionStageAdmissionCertificate] = {}
    for row in cert_rows:
        if row.stage in cert_by_stage:
            raise ValueError(f"duplicate admitted Section 9 stage {row.stage}")
        if not row.formal_admission_ready:
            raise ValueError(f"stage {row.stage} is not a complete formal admission")
        if row.q_big != hypothesis.q_big:
            raise ValueError("all admitted stages must use the hypothesis common q_big")
        expected_lower = hypothesis.first_scale * (1 << (row.stage - 1))
        if row.schedule_lower_bound != expected_lower:
            raise ValueError("stage certificate was derived from a different schedule first_scale")
        if row.scale < expected_lower:
            raise ArithmeticError("admitted stage scale violates the supplied schedule lower bound")
        cert_by_stage[row.stage] = row

    stages = tuple(sorted(cert_by_stage))
    expected_stages = tuple(range(1, len(stages) + 1))
    if stages != expected_stages:
        raise ValueError(
            "finite Eq. (9.21) truncation must be a contiguous positive-stage prefix 1..N"
        )

    values_by_stage: dict[int, Section9CorrectionStageValue] = {}
    for row in value_rows:
        if row.stage in values_by_stage:
            raise ValueError(f"duplicate Section 9 stage value {row.stage}")
        values_by_stage[row.stage] = row
    if set(values_by_stage) != set(stages):
        missing = sorted(set(stages) - set(values_by_stage))
        extra = sorted(set(values_by_stage) - set(stages))
        raise ValueError(
            "stage values must cover exactly the admitted prefix; "
            f"missing={missing}, extra={extra}"
        )

    A_sum = np.zeros(3, dtype=float)
    B_sum = 0.0
    p_sum = 0.0
    contributions: list[Section9FinitePrefixStageContribution] = []

    for stage in stages:
        cert = cert_by_stage[stage]
        value = values_by_stage[stage]
        argument = cert.scale * q_value
        short_circuit = argument >= 1
        if short_circuit:
            weight = 0.0
        else:
            weight = _finite_scalar(
                cutoff_evaluator(float(argument)),
                f"cutoff_evaluator(a_{stage} q)",
            )
            if not 0.0 <= weight <= 1.0:
                raise ValueError("cutoff_evaluator must return a value in [0,1]")

        A_term = np.array(weight * value.A, dtype=float, copy=True)
        A_term.setflags(write=False)
        B_term = weight * value.B
        p_term = weight * value.p
        A_sum += A_term
        B_sum += B_term
        p_sum += p_term
        contributions.append(
            Section9FinitePrefixStageContribution(
                stage=stage,
                scale=cert.scale,
                cutoff_argument=argument,
                cutoff_weight=weight,
                A_term=A_term,
                B_term=B_term,
                p_term=p_term,
                value_provenance=value.provenance,
                support_zero_short_circuit=short_circuit,
            )
        )

    if not (np.all(np.isfinite(A_sum)) and math.isfinite(B_sum) and math.isfinite(p_sum)):
        raise ArithmeticError("finite-prefix accumulation overflowed")
    A_sum.setflags(write=False)
    result = Section9FinitePrefixEvaluationCertificate(
        q=q_value,
        q_big=hypothesis.q_big,
        prefix_order=len(stages),
        stages=stages,
        contributions=tuple(contributions),
        A_correction_prefix=A_sum,
        B_correction_prefix=float(B_sum),
        p_correction_prefix=float(p_sum),
        cutoff_evaluator_provenance=cutoff_provenance,
    )
    if not result.formal_prefix_ready:
        raise ArithmeticError("Section 9 finite-prefix truth-status invariant failed")
    return result
