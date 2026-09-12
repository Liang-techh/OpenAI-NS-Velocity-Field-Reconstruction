"""Bind the actual post-temporal debt to the compact Section 8 rank repair.

The pinned Lean development proves
``ActualIntermediateDebtBounds.afterTemporal_debt_from_stepData``: checked
correction-step data imply an ``UnweightedClass`` bound for the full
three-component ``CorrectionState.debt`` of the literal ``postTemporal`` state.
That theorem is the formal handoff immediately before the actual compact rank
repair.

This module does not materialize that debt.  Instead it admits a content-free,
fail-closed ``lean-formal-export`` witness and composes it with the already
landed ``CompactMeanCorrectionAdmission``.  The composition requires the exact
same context, post-temporal state, band/slow point, and debt instance.

The finite-head signed cross defect is intentionally *not* identified with this
state debt.  That equality needs its own formal theorem/export; sampled or
inferred data are never accepted here.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .compact_mean_correction_formal import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
    CompactMeanCorrectionAdmission,
)


PINNED_INTERMEDIATE_DEBT_FILE = "NavierStokes/ActualIntermediateDebtBounds.lean"
PINNED_CORRECTION_ANALYTIC_STEP_FILE = "NavierStokes/CorrectionAnalyticStep.lean"
PINNED_CORRECTION_STATE_FILE = "NavierStokes/CorrectionState.lean"
PINNED_AFTER_TEMPORAL_DEBT_THEOREM = (
    "NavierStokes.ActualIntermediateDebtBounds.afterTemporal_debt_from_stepData"
)
PINNED_STAGE_DEBT_THEOREM = (
    "NavierStokes.ActualIntermediateDebtBounds.stage_debt_from_stepData"
)
PINNED_DEBT_SYMBOL = "NavierStokes.CorrectionState.debt"
REQUIRED_FORMAL_FILES = (
    PINNED_INTERMEDIATE_DEBT_FILE,
    PINNED_CORRECTION_ANALYTIC_STEP_FILE,
    PINNED_CORRECTION_STATE_FILE,
)
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_AFTER_TEMPORAL_DEBT_THEOREM,
    PINNED_STAGE_DEBT_THEOREM,
    PINNED_DEBT_SYMBOL,
)
PRODUCER_KIND = "lean-formal-export"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


@dataclass(frozen=True)
class ActualIntermediateRankDebtWitness:
    """Opaque export metadata for one post-temporal debt family and evaluation.

    ``debt_family_id`` names the exact Lean function
    ``CorrectionState.debt (ActualPrimary.commonContext B) (postTemporal x)``.
    ``band_debt_id`` names its evaluation at ``band_index`` and
    ``slow_point_id``.  Python never reconstructs either object.
    """

    B: int
    N0: int
    cycle_state_id: str
    static_data_id: str
    cycle_invariant_id: str
    particular_inputs_id: str
    step_data_id: str
    sigma_repr: str
    context_id: str
    post_temporal_state_id: str
    debt_family_id: str
    band_index: int
    slow_point_id: str
    band_debt_id: str
    application_id: str
    producer_kind: str
    provenance: str

    same_B_N0_certified: bool
    same_cycle_state_certified: bool
    common_context_identity_certified: bool
    post_temporal_state_identity_certified: bool
    step_data_for_invariant_certified: bool
    sigma_lower_bound_certified: bool
    debt_family_identity_certified: bool
    band_debt_evaluation_identity_certified: bool
    theorem_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_INTERMEDIATE_DEBT_FILE
    theorem_symbol: str = PINNED_AFTER_TEMPORAL_DEBT_THEOREM
    dependency_files: tuple[str, ...] = REQUIRED_FORMAL_FILES
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "band_index", _natural(self.band_index, "band_index"))

        for name in (
            "cycle_state_id",
            "static_data_id",
            "cycle_invariant_id",
            "particular_inputs_id",
            "step_data_id",
            "sigma_repr",
            "context_id",
            "post_temporal_state_id",
            "debt_family_id",
            "slow_point_id",
            "band_debt_id",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_INTERMEDIATE_DEBT_FILE:
            raise ValueError("formal_file must match ActualIntermediateDebtBounds.lean exactly")
        if self.theorem_symbol != PINNED_AFTER_TEMPORAL_DEBT_THEOREM:
            raise ValueError("theorem_symbol must be afterTemporal_debt_from_stepData")
        if self.dependency_files != REQUIRED_FORMAL_FILES:
            raise ValueError("dependency_files must match the pinned intermediate-debt chain exactly")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned intermediate-debt chain exactly")

        for name in (
            "same_B_N0_certified",
            "same_cycle_state_certified",
            "common_context_identity_certified",
            "post_temporal_state_identity_certified",
            "step_data_for_invariant_certified",
            "sigma_lower_bound_certified",
            "debt_family_identity_certified",
            "band_debt_evaluation_identity_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class ActualIntermediateRankDebtAdmission:
    """Compose the actual post-temporal debt export with compact rank repair."""

    intermediate: ActualIntermediateRankDebtWitness
    compact_rank: CompactMeanCorrectionAdmission

    def __post_init__(self) -> None:
        if not isinstance(self.intermediate, ActualIntermediateRankDebtWitness):
            raise TypeError("intermediate must be an ActualIntermediateRankDebtWitness")
        if not isinstance(self.compact_rank, CompactMeanCorrectionAdmission):
            raise TypeError("compact_rank must be a CompactMeanCorrectionAdmission")

        rank = self.compact_rank.witness
        if self.intermediate.formal_repository != rank.formal_repository:
            raise ValueError("intermediate debt and compact rank must use the same formal repository")
        if self.intermediate.formal_commit != rank.formal_commit:
            raise ValueError("intermediate debt and compact rank must use the same formal commit")
        if self.intermediate.context_id != rank.context_id:
            raise ValueError("compact rank context is not the actual post-temporal context")
        if self.intermediate.post_temporal_state_id != rank.state_id:
            raise ValueError("compact rank state is not the actual post-temporal state")
        if self.intermediate.band_index != rank.band_index:
            raise ValueError("compact rank band does not match the admitted debt evaluation")
        if self.intermediate.slow_point_id != rank.slow_point_id:
            raise ValueError("compact rank slow point does not match the admitted debt evaluation")
        if self.intermediate.band_debt_id != rank.debt_id:
            raise ValueError("compact rank debt is not the admitted post-temporal debt evaluation")

    @property
    def after_temporal_debt_theorem_admitted(self) -> bool:
        return True

    @property
    def actual_post_temporal_debt_bound_to_compact_rank(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def theorem_application_machine_replayed(self) -> bool:
        return False

    @property
    def actual_cycle_state_materialized(self) -> bool:
        return False

    @property
    def actual_post_temporal_debt_values_materialized(self) -> bool:
        return False

    @property
    def finite_head_signed_defect_link_verified(self) -> bool:
        return False

    @property
    def compact_mean_correction_field_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
