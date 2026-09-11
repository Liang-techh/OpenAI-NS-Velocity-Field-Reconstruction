"""Fail-closed bridge from a physical Section 10 slab to a complete Section 9 prefix.

The landed Section 9 machinery now contains three separate facts:

* an infinite doubling-schedule/local-finiteness theorem input,
* an exact Eq. (4.1) enclosure of a closed pre-endpoint Section 10 slab by a
  compact q-strip, and
* a numerical evaluator for a contiguous finite prefix of admitted Eq. (9.21)
  correction stages.

This module composes those facts without upgrading their provenance.  If the
finite prefix reaches the slab certificate's ``max_potentially_active_stage``,
then every omitted positive stage is forced into the paper cutoff's support-zero
region on the *whole slab*.  Thus the prefix is structurally stage-complete on
that slab.  The statement concerns support only: caller-supplied correction
values and a generic pointwise cutoff evaluator remain unverified, the endpoint
``t=1`` is not covered, and the infinite Eq. (9.21) field is not thereby
constructed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_finite_prefix import Section9FinitePrefixEvaluationCertificate
from .section9_local_sum_admission import Section9LocalSumHypothesis
from .section9_physical_window import Section9PhysicalSlabLocalSumCertificate


@dataclass(frozen=True)
class Section9PhysicalPrefixCompletionCertificate:
    """Support-completeness of one evaluated prefix on one pre-endpoint slab."""

    t_lower: Fraction
    t_upper: Fraction
    q_lower_bound: Fraction
    q_upper_bound: Fraction
    evaluation_q: Fraction
    prefix_order: int
    max_potentially_active_stage: int
    first_guaranteed_inactive_stage: int
    first_omitted_stage: int
    omitted_tail_support_zero_on_whole_slab: bool = True
    evaluated_q_inside_slab_enclosure: bool = True
    admitted_prefix_matches_evaluation: bool = True
    physical_schedule_coherence_verified: bool = True
    endpoint_covered: bool = False
    actual_correction_field_values_verified: bool = False
    paper_fixed_cutoff_pointwise_evaluator_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    endpoint_uniform_residual_majorants_verified: bool = False
    proposition_9_9_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_slab_prefix_ready(self) -> bool:
        return (
            self.prefix_order >= self.max_potentially_active_stage
            and self.first_omitted_stage == self.prefix_order + 1
            and self.first_guaranteed_inactive_stage
            == self.max_potentially_active_stage + 1
            and self.omitted_tail_support_zero_on_whole_slab
            and self.evaluated_q_inside_slab_enclosure
            and self.admitted_prefix_matches_evaluation
            and self.physical_schedule_coherence_verified
            and not self.endpoint_covered
            and not self.actual_correction_field_values_verified
            and not self.paper_fixed_cutoff_pointwise_evaluator_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.endpoint_uniform_residual_majorants_verified
            and not self.proposition_9_9_verified
            and not self.paper_exact_velocity_available
        )


def certify_section9_physical_prefix_completion(
    hypothesis: Section9LocalSumHypothesis,
    physical_slab: Section9PhysicalSlabLocalSumCertificate,
    evaluation: Section9FinitePrefixEvaluationCertificate,
    stage_certificates: Iterable[Section9CorrectionStageAdmissionCertificate],
) -> Section9PhysicalPrefixCompletionCertificate:
    """Certify that no omitted positive stage can survive on ``physical_slab``.

    ``evaluation`` is still only a pointwise finite-prefix evaluation.  This
    function does not assert that its numerical values are the paper corrections.
    It only checks that the evaluated contiguous prefix contains every stage that
    the independently certified infinite doubling schedule allows to be nonzero
    anywhere on the whole physical slab.
    """
    if not isinstance(hypothesis, Section9LocalSumHypothesis):
        raise TypeError("hypothesis must be a Section9LocalSumHypothesis")
    if not isinstance(physical_slab, Section9PhysicalSlabLocalSumCertificate):
        raise TypeError("physical_slab must be a Section9PhysicalSlabLocalSumCertificate")
    if not isinstance(evaluation, Section9FinitePrefixEvaluationCertificate):
        raise TypeError("evaluation must be a Section9FinitePrefixEvaluationCertificate")
    if not evaluation.formal_prefix_ready:
        raise ValueError("evaluation must be a complete formal finite-prefix evaluation")

    local_sum = physical_slab.local_sum
    if (
        local_sum.q_big != hypothesis.q_big
        or local_sum.first_scale != hypothesis.first_scale
        or local_sum.evidence_kind != hypothesis.kind
        or local_sum.evidence_provenance != hypothesis.provenance
    ):
        raise ValueError("physical slab certificate must come from the supplied local-sum hypothesis")
    if evaluation.q_big != hypothesis.q_big:
        raise ValueError("finite-prefix evaluation must use the supplied common q_big")
    if not physical_slab.q_lower_bound <= evaluation.q <= physical_slab.q_upper_bound:
        raise ValueError("finite-prefix evaluation q lies outside the physical slab q enclosure")

    rows = tuple(stage_certificates)
    if not rows:
        raise ValueError("stage certificates must cover the evaluated prefix")
    if not all(isinstance(row, Section9CorrectionStageAdmissionCertificate) for row in rows):
        raise TypeError("every stage certificate must be a Section9CorrectionStageAdmissionCertificate")

    by_stage: dict[int, Section9CorrectionStageAdmissionCertificate] = {}
    for row in rows:
        if row.stage in by_stage:
            raise ValueError(f"duplicate admitted Section 9 stage {row.stage}")
        if not row.formal_admission_ready:
            raise ValueError(f"stage {row.stage} is not a complete formal admission")
        by_stage[row.stage] = row
    if tuple(sorted(by_stage)) != evaluation.stages:
        raise ValueError("stage certificates must cover exactly the evaluated contiguous prefix")

    for contribution in evaluation.contributions:
        row = by_stage[contribution.stage]
        expected_lower = hypothesis.first_scale * (1 << (row.stage - 1))
        if row.q_big != hypothesis.q_big:
            raise ValueError("admitted stage uses a different common q_big")
        if row.schedule_lower_bound != expected_lower:
            raise ValueError("admitted stage was derived from a different schedule first_scale")
        if row.scale != contribution.scale:
            raise ValueError("evaluated stage scale does not match its admission certificate")

    required = physical_slab.max_potentially_active_stage
    if evaluation.prefix_order < required:
        raise ValueError(
            "evaluated prefix is too short for whole-slab support completeness: "
            f"need stages 1..{required}"
        )

    result = Section9PhysicalPrefixCompletionCertificate(
        t_lower=physical_slab.t_lower,
        t_upper=physical_slab.t_upper,
        q_lower_bound=physical_slab.q_lower_bound,
        q_upper_bound=physical_slab.q_upper_bound,
        evaluation_q=evaluation.q,
        prefix_order=evaluation.prefix_order,
        max_potentially_active_stage=required,
        first_guaranteed_inactive_stage=physical_slab.first_guaranteed_inactive_stage,
        first_omitted_stage=evaluation.prefix_order + 1,
    )
    if not result.formal_slab_prefix_ready:
        raise ArithmeticError("Section 9 physical-prefix completion invariant failed")
    return result
