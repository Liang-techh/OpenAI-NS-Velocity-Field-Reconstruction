"""Fail-closed Section 8 -> Section 9 rank-input identity bridge.

The pinned formal chain uses exactly the same post-temporal debt twice:

* ``LocalRankDefect.RankGeometry.solved_rows`` cancels that debt by the actual
  rank increment and ``debt_eq_remainders`` identifies the post-rank debt; and
* ``ActualStageEstimates.RunData.rank_class`` obtains the native rank-class
  estimate from ``afterTemporal_debt_from_stepData`` and feeds it to
  ``ActualStageEstimates.rankInput`` / ``actualCycleRankInput``.

This module binds those opaque formal objects to one another and to an already
admitted ``Section9ActualStageEstimatesAdmission``.  It deliberately does not
materialize the noncomputable rank increment, the physical rank family, or any
Section 9 field/value.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .section8_local_rank_solve import Section8LocalRankSolveAdmission
from .section9_actual_stage_estimates import Section9ActualStageEstimatesAdmission


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_STAGE_FILE = "NavierStokes/ActualStageEstimates.lean"
PINNED_PHYSICAL_FILE = "NavierStokes/ActualPhysicalStageBounds.lean"
PINNED_INTERMEDIATE_FILE = "NavierStokes/ActualIntermediateDebtBounds.lean"
PINNED_LOCAL_RANK_FILE = "NavierStokes/LocalRankDefect.lean"
PINNED_RANK_CLASS = "NavierStokes.ActualStageEstimates.RunData.rank_class"
PINNED_RANK_INPUT = "NavierStokes.ActualStageEstimates.rankInput"
PINNED_ACTUAL_CYCLE_RANK_INPUT = (
    "NavierStokes.ActualPhysicalStageBounds.actualCycleRankInput"
)
PINNED_AFTER_TEMPORAL_DEBT = (
    "NavierStokes.ActualIntermediateDebtBounds.afterTemporal_debt_from_stepData"
)
PINNED_SOLVED_ROWS = "NavierStokes.LocalRankDefect.RankGeometry.solved_rows"
PINNED_DEBT_REMAINDERS = (
    "NavierStokes.LocalRankDefect.RankGeometry.debt_eq_remainders"
)
REQUIRED_FORMAL_FILES = (
    PINNED_STAGE_FILE,
    PINNED_PHYSICAL_FILE,
    PINNED_INTERMEDIATE_FILE,
    PINNED_LOCAL_RANK_FILE,
)
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_RANK_CLASS,
    PINNED_RANK_INPUT,
    PINNED_ACTUAL_CYCLE_RANK_INPUT,
    PINNED_AFTER_TEMPORAL_DEBT,
    PINNED_SOLVED_ROWS,
    PINNED_DEBT_REMAINDERS,
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
class Section8Section9RankHandoffWitness:
    """Opaque identities for one exact rank-input handoff.

    ``cycle_state_id`` / ``post_temporal_state_id`` / ``debt_family_id`` bind
    to the already-landed Section 8 actual intermediate-debt chain.
    ``rank_stage_state_id`` and ``remainder_family_id`` bind to its local rank
    solve.  ``run_data_id`` and ``stage_estimates_id`` bind to the Section 9
    actual finite-stage admission.  Python never reconstructs any of them.
    """

    B: int
    N0: int
    stage_index: int
    cycle_state_id: str
    post_temporal_state_id: str
    debt_family_id: str
    rank_stage_state_id: str
    remainder_family_id: str
    run_data_id: str
    rank_class_id: str
    rank_input_id: str
    stage_estimates_id: str
    rank_class_application_id: str
    rank_input_application_id: str
    producer_kind: str
    provenance: str

    same_cycle_state_certified: bool
    same_post_temporal_state_certified: bool
    same_debt_family_certified: bool
    rank_class_uses_after_temporal_debt_certified: bool
    rank_input_uses_rank_class_certified: bool
    local_rank_solve_uses_same_debt_certified: bool
    post_rank_remainder_identity_certified: bool
    same_run_data_certified: bool
    same_stage_estimates_certified: bool
    handoff_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_files: tuple[str, ...] = REQUIRED_FORMAL_FILES
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "stage_index", _natural(self.stage_index, "stage_index"))
        for name in (
            "cycle_state_id",
            "post_temporal_state_id",
            "debt_family_id",
            "rank_stage_state_id",
            "remainder_family_id",
            "run_data_id",
            "rank_class_id",
            "rank_input_id",
            "stage_estimates_id",
            "rank_class_application_id",
            "rank_input_application_id",
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
        if self.formal_files != REQUIRED_FORMAL_FILES:
            raise ValueError("formal_files must match the pinned rank-input chain exactly")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned rank-input chain exactly")
        for name in (
            "same_cycle_state_certified",
            "same_post_temporal_state_certified",
            "same_debt_family_certified",
            "rank_class_uses_after_temporal_debt_certified",
            "rank_input_uses_rank_class_certified",
            "local_rank_solve_uses_same_debt_certified",
            "post_rank_remainder_identity_certified",
            "same_run_data_certified",
            "same_stage_estimates_certified",
            "handoff_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section8Section9RankHandoffAdmission:
    """Compose the landed Section 8 local rank solve with Section 9 rank input."""

    local_rank: Section8LocalRankSolveAdmission
    stage_estimates: Section9ActualStageEstimatesAdmission
    witness: Section8Section9RankHandoffWitness

    def __post_init__(self) -> None:
        if not isinstance(self.local_rank, Section8LocalRankSolveAdmission):
            raise TypeError("local_rank must be a Section8LocalRankSolveAdmission")
        if not isinstance(self.stage_estimates, Section9ActualStageEstimatesAdmission):
            raise TypeError(
                "stage_estimates must be a Section9ActualStageEstimatesAdmission"
            )
        if not isinstance(self.witness, Section8Section9RankHandoffWitness):
            raise TypeError("witness must be a Section8Section9RankHandoffWitness")

        intermediate = self.local_rank.intermediate_rank.intermediate
        rank = self.local_rank.witness
        stage = self.stage_estimates.witness

        if (self.witness.B, self.witness.N0) != (intermediate.B, intermediate.N0):
            raise ValueError("rank handoff B/N0 does not match the actual correction cycle")
        if self.witness.formal_repository != intermediate.formal_repository:
            raise ValueError("rank handoff repository does not match Section 8")
        if self.witness.formal_commit != intermediate.formal_commit:
            raise ValueError("rank handoff commit does not match Section 8")
        if self.witness.formal_repository != stage.formal_repository:
            raise ValueError("rank handoff repository does not match Section 9")
        if self.witness.formal_commit != stage.formal_commit:
            raise ValueError("rank handoff commit does not match Section 9")
        if self.witness.cycle_state_id != intermediate.cycle_state_id:
            raise ValueError("rank handoff cycle state is cross-wired")
        if self.witness.post_temporal_state_id != intermediate.post_temporal_state_id:
            raise ValueError("rank handoff post-temporal state is cross-wired")
        if self.witness.debt_family_id != intermediate.debt_family_id:
            raise ValueError("rank handoff debt family is cross-wired")
        if self.witness.rank_stage_state_id != rank.rank_stage_state_id:
            raise ValueError("rank handoff rank-stage state is cross-wired")
        if self.witness.remainder_family_id != rank.remainder_family_id:
            raise ValueError("rank handoff remainder family is cross-wired")
        if self.witness.run_data_id != stage.run_data_id:
            raise ValueError("rank handoff RunData is not the admitted Section 9 RunData")
        if self.witness.stage_estimates_id != stage.stage_estimates_id:
            raise ValueError("rank handoff StageEstimates object is cross-wired")

    @property
    def same_after_temporal_debt_bound_to_section8_and_section9(self) -> bool:
        return True

    @property
    def section9_rank_class_admitted(self) -> bool:
        return True

    @property
    def section9_rank_input_identity_admitted(self) -> bool:
        return True

    @property
    def post_rank_remainder_identity_preserved(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def theorem_chain_machine_replayed(self) -> bool:
        return False

    @property
    def actual_rank_increment_materialized(self) -> bool:
        return False

    @property
    def actual_rank_stage_state_materialized(self) -> bool:
        return False

    @property
    def actual_section9_rank_input_materialized(self) -> bool:
        return False

    @property
    def section9_iteration_machine_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
