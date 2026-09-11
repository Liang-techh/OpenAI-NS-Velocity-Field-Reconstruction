"""Fail-closed admission for smooth Section 9 correction extensions.

Proposition 9.9 / Eq. (9.21) ultimately sums the positive-stage corrections
through the fixed factors ``chi(a_j q)``.  The already-landed local-sum module
checks only the *arithmetic* local-finiteness consequence of an externally
certified common q-domain and infinite doubling schedule.  It does not certify
that a concrete correction at stage ``j`` is one of the paper's smooth
zero-extended fields.

This module supplies that missing structural gate.  For one positive stage it
requires separate theorem/proof witnesses for the three Eq. (9.21) components
``A_j``, ``B_j`` and ``p_j``.  Every witness must certify the same common
q-domain, the same cutoff scale ``a_j``, membership in the admitted schedule,
the paper cutoff-support rule, and smooth zero extension.  The module then
checks component completeness and the exact schedule lower bound

    a_j >= a_1 * 2**(j-1).

It does not inspect or fabricate field values, infer schedule membership from
that lower bound, construct the Eq. (9.21) sum, or prove Proposition 9.9.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Rational, Real
import math
from typing import Iterable

from .section9_local_sum_admission import Section9LocalSumHypothesis


_REQUIRED_COMPONENTS = ("A", "B", "p")
_RIGOROUS_KINDS = frozenset({
    "paper-derived",
    "lean-derived",
    "certified-numerical",
    "rigorous-external",
})


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
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


def _require_true(value: object, name: str) -> None:
    if value is not True:
        raise ValueError(f"{name} must be an independently justified theorem/proof fact")


@dataclass(frozen=True)
class Section9CorrectionExtensionWitness:
    """External theorem/proof witness for one Eq. (9.21) correction component.

    ``component`` is exactly one of ``A``, ``B`` or ``p``.  The four boolean
    fields are deliberately assertions supplied by an independent theorem or
    proof artifact; sampled field values are not accepted as substitutes.
    ``schedule_membership_verified`` means the supplied ``scale`` is genuinely
    the paper's ``a_j`` for this stage, not merely an arbitrary integer that
    happens to satisfy the numerical lower bound checked later.
    """

    stage: int
    component: str
    q_big: Fraction | float | int
    scale: int
    kind: str
    provenance: str
    common_domain_verified: bool
    schedule_membership_verified: bool
    cutoff_support_verified: bool
    smooth_zero_extension_verified: bool

    def __post_init__(self) -> None:
        stage = _positive_int(self.stage, "stage")
        scale = _positive_int(self.scale, "scale")
        q_big = _fraction(self.q_big, "q_big")
        component = str(self.component)
        if component not in _REQUIRED_COMPONENTS:
            raise ValueError(f"component must be one of {_REQUIRED_COMPONENTS}")
        if q_big > 1:
            raise ValueError("Section 9 common q-domain requires q_big<=1")
        if self.kind not in _RIGOROUS_KINDS:
            raise ValueError(f"kind must be one of {sorted(_RIGOROUS_KINDS)}")
        provenance = str(self.provenance).strip()
        if not provenance:
            raise ValueError("provenance must be nonempty")
        _require_true(self.common_domain_verified, "common_domain_verified")
        _require_true(self.schedule_membership_verified, "schedule_membership_verified")
        _require_true(self.cutoff_support_verified, "cutoff_support_verified")
        _require_true(self.smooth_zero_extension_verified, "smooth_zero_extension_verified")
        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "scale", scale)
        object.__setattr__(self, "q_big", q_big)
        object.__setattr__(self, "component", component)
        object.__setattr__(self, "provenance", provenance)


@dataclass(frozen=True)
class Section9CorrectionStageAdmissionCertificate:
    """Coherent formal admission of all three correction types at one stage."""

    stage: int
    q_big: Fraction
    scale: int
    schedule_lower_bound: int
    components: tuple[str, ...]
    evidence_kinds: tuple[str, ...]
    evidence_provenance: tuple[str, ...]
    common_domain_verified: bool = True
    schedule_membership_witnessed: bool = True
    schedule_lower_bound_verified: bool = True
    cutoff_support_verified: bool = True
    smooth_zero_extension_verified: bool = True
    component_completeness_verified: bool = True
    source_theorems_machine_verified: bool = False
    actual_correction_field_values_verified: bool = False
    infinite_schedule_constructed_here: bool = False
    eq_9_21_sum_constructed: bool = False
    proposition_9_9_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_admission_ready(self) -> bool:
        return (
            self.components == _REQUIRED_COMPONENTS
            and self.scale >= self.schedule_lower_bound
            and self.common_domain_verified
            and self.schedule_membership_witnessed
            and self.schedule_lower_bound_verified
            and self.cutoff_support_verified
            and self.smooth_zero_extension_verified
            and self.component_completeness_verified
            and not self.source_theorems_machine_verified
            and not self.actual_correction_field_values_verified
            and not self.infinite_schedule_constructed_here
            and not self.eq_9_21_sum_constructed
            and not self.proposition_9_9_verified
            and not self.paper_exact_velocity_available
        )


def admit_section9_correction_stage(
    hypothesis: Section9LocalSumHypothesis,
    witnesses: Iterable[Section9CorrectionExtensionWitness],
) -> Section9CorrectionStageAdmissionCertificate:
    """Admit one positive Eq. (9.21) stage after theorem-level extension checks.

    The global ``hypothesis`` remains an external theorem input.  This function
    checks that the three component witnesses are coherent with it and verifies
    the exact consequence ``a_j >= a_1 2**(j-1)``.  It intentionally does not
    treat that consequence as proof that an arbitrary scale belongs to the
    schedule; each witness must separately certify schedule membership.
    """

    if not isinstance(hypothesis, Section9LocalSumHypothesis):
        raise TypeError("hypothesis must be a Section9LocalSumHypothesis")
    rows = tuple(witnesses)
    if not rows:
        raise ValueError("correction-stage witnesses must be nonempty")
    if not all(isinstance(row, Section9CorrectionExtensionWitness) for row in rows):
        raise TypeError("every entry must be a Section9CorrectionExtensionWitness")

    by_component: dict[str, Section9CorrectionExtensionWitness] = {}
    for row in rows:
        if row.component in by_component:
            raise ValueError(f"duplicate Section 9 component {row.component!r}")
        by_component[row.component] = row
    supplied = set(by_component)
    required = set(_REQUIRED_COMPONENTS)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            "one Section 9 stage must cover exactly A, B and p; "
            f"missing={missing}, extra={extra}"
        )

    ordered = tuple(by_component[name] for name in _REQUIRED_COMPONENTS)
    first = ordered[0]
    for row in ordered[1:]:
        if row.stage != first.stage:
            raise ValueError("all correction components must use one stage")
        if row.q_big != first.q_big:
            raise ValueError("all correction components must use one common q_big")
        if row.scale != first.scale:
            raise ValueError("all correction components must use the same stage cutoff scale")
    if first.q_big != hypothesis.q_big:
        raise ValueError("correction witnesses must use the local-sum hypothesis common q_big")

    lower_bound = hypothesis.first_scale * (1 << (first.stage - 1))
    if first.scale < lower_bound:
        raise ValueError(
            "witnessed stage scale violates the admitted doubling-schedule lower bound: "
            f"a_{first.stage}={first.scale} < {lower_bound}"
        )

    certificate = Section9CorrectionStageAdmissionCertificate(
        stage=first.stage,
        q_big=first.q_big,
        scale=first.scale,
        schedule_lower_bound=lower_bound,
        components=_REQUIRED_COMPONENTS,
        evidence_kinds=tuple(row.kind for row in ordered),
        evidence_provenance=tuple(row.provenance for row in ordered),
    )
    if not certificate.formal_admission_ready:
        raise ArithmeticError("Section 9 correction-extension admission invariant failed")
    return certificate
