"""Bind the actual post-temporal rank debt to the local Section 8 solve.

The pinned formal development goes beyond the isolated five-row algebra already
admitted by :mod:`compact_mean_correction_formal`.  In
``NavierStokes/LocalRankDefect.lean``, ``RankGeometry.solved_rows`` proves on the
actual open slow domain that the variable-gauge rank increment satisfies

``linearRows ... rankIncrementState = -CorrectionState.debt ...``.

The same formal chain proves preservation of the two required masses,
``RankGeometry.debt_eq_remainders`` for the updated rank-stage state, and common
moving-support retention.  These are the theorem-facing handoff from the
Section 8 compact rank correction to the Section 9 residual iteration.

This module is deliberately content-free.  It accepts only exact
``lean-formal-export`` identities and composes them with the already-landed
``ActualIntermediateRankDebtAdmission``.  It does not reconstruct the opaque
Lean bump, materialize the rank increment, infer debt values, or claim that a
finite numerical solve is the paper construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .actual_intermediate_rank_debt import ActualIntermediateRankDebtAdmission
from .compact_mean_correction_formal import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
)


PINNED_LOCAL_RANK_FILE = "NavierStokes/LocalRankDefect.lean"
PINNED_CORRECTION_STATE_FILE = "NavierStokes/CorrectionState.lean"
PINNED_MEAN_RANK_FILE = "NavierStokes/MeanRankUpdate.lean"
PINNED_SOLVED_ROWS_THEOREM = "NavierStokes.LocalRankDefect.RankGeometry.solved_rows"
PINNED_PRESERVE_MASSES_THEOREM = (
    "NavierStokes.LocalRankDefect.RankGeometry.preserve_masses"
)
PINNED_DEBT_REMAINDERS_THEOREM = (
    "NavierStokes.LocalRankDefect.RankGeometry.debt_eq_remainders"
)
PINNED_SUPPORT_THEOREM = (
    "NavierStokes.LocalRankDefect.RankGeometry.increment_supportedGauge"
)
PINNED_FIVE_ROWS_THEOREM = "NavierStokes.LocalRankDefect.RankGeometry.fiveRows"
PINNED_RANK_GEOMETRY = "NavierStokes.LocalRankDefect.RankGeometry"
REQUIRED_FORMAL_FILES = (
    PINNED_LOCAL_RANK_FILE,
    PINNED_CORRECTION_STATE_FILE,
    PINNED_MEAN_RANK_FILE,
)
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_RANK_GEOMETRY,
    PINNED_FIVE_ROWS_THEOREM,
    PINNED_SOLVED_ROWS_THEOREM,
    PINNED_PRESERVE_MASSES_THEOREM,
    PINNED_DEBT_REMAINDERS_THEOREM,
    PINNED_SUPPORT_THEOREM,
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
class Section8LocalRankSolveWitness:
    """Opaque formal-export metadata for one actual local rank-stage solve.

    ``input_state_id`` is the literal ``postTemporal x`` state already bound by
    :class:`ActualIntermediateRankDebtAdmission`.  ``increment_state_id`` names
    ``VariableGaugeMean.rankIncrementState`` at that state and
    ``rank_stage_state_id`` names the corresponding ``rankStageState``.
    ``remainder_family_id`` names the exact ``LocalRankDefect.remainders``
    expression to which the post-rank debt is reduced.
    """

    B: int
    N0: int
    rank_data_id: str
    gauge_data_id: str
    local_domain_id: str
    context_id: str
    input_state_id: str
    debt_family_id: str
    band_index: int
    slow_point_id: str
    axial_direction_id: str
    increment_state_id: str
    rank_stage_state_id: str
    remainder_family_id: str
    application_id: str
    producer_kind: str
    provenance: str

    rank_geometry_certified: bool
    local_domain_open_certified: bool
    support_window_certified: bool
    local_operators_certified: bool
    base_smooth_certified: bool
    angular_base_slow_certified: bool
    axial_base_slow_certified: bool
    current_mean_local_triple_certified: bool
    covariance_local_shell_certified: bool
    five_rows_application_certified: bool
    solved_rows_application_certified: bool
    preserve_masses_application_certified: bool
    debt_eq_remainders_application_certified: bool
    support_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_LOCAL_RANK_FILE
    theorem_symbol: str = PINNED_SOLVED_ROWS_THEOREM
    dependency_files: tuple[str, ...] = REQUIRED_FORMAL_FILES
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "band_index", _natural(self.band_index, "band_index"))

        for name in (
            "rank_data_id",
            "gauge_data_id",
            "local_domain_id",
            "context_id",
            "input_state_id",
            "debt_family_id",
            "slow_point_id",
            "axial_direction_id",
            "increment_state_id",
            "rank_stage_state_id",
            "remainder_family_id",
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
        if self.formal_file != PINNED_LOCAL_RANK_FILE:
            raise ValueError("formal_file must match LocalRankDefect.lean exactly")
        if self.theorem_symbol != PINNED_SOLVED_ROWS_THEOREM:
            raise ValueError("theorem_symbol must be RankGeometry.solved_rows")
        if self.dependency_files != REQUIRED_FORMAL_FILES:
            raise ValueError("dependency_files must match the pinned local-rank chain exactly")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned local-rank chain exactly")

        for name in (
            "rank_geometry_certified",
            "local_domain_open_certified",
            "support_window_certified",
            "local_operators_certified",
            "base_smooth_certified",
            "angular_base_slow_certified",
            "axial_base_slow_certified",
            "current_mean_local_triple_certified",
            "covariance_local_shell_certified",
            "five_rows_application_certified",
            "solved_rows_application_certified",
            "preserve_masses_application_certified",
            "debt_eq_remainders_application_certified",
            "support_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section8LocalRankSolveAdmission:
    """Compose the actual post-temporal debt with the local five-row solve."""

    intermediate_rank: ActualIntermediateRankDebtAdmission
    witness: Section8LocalRankSolveWitness

    def __post_init__(self) -> None:
        if not isinstance(self.intermediate_rank, ActualIntermediateRankDebtAdmission):
            raise TypeError("intermediate_rank must be an ActualIntermediateRankDebtAdmission")
        if not isinstance(self.witness, Section8LocalRankSolveWitness):
            raise TypeError("witness must be a Section8LocalRankSolveWitness")

        intermediate = self.intermediate_rank.intermediate
        rank = self.intermediate_rank.compact_rank.witness

        if (self.witness.B, self.witness.N0) != (intermediate.B, intermediate.N0):
            raise ValueError("local rank solve B/N0 does not match the actual correction cycle")
        if self.witness.formal_repository != intermediate.formal_repository:
            raise ValueError("local rank solve repository does not match the intermediate debt")
        if self.witness.formal_commit != intermediate.formal_commit:
            raise ValueError("local rank solve commit does not match the intermediate debt")
        if self.witness.rank_data_id != rank.rank_data_id:
            raise ValueError("local rank solve rank data does not match compact rank repair")
        if self.witness.context_id != intermediate.context_id:
            raise ValueError("local rank solve context is not the actual post-temporal context")
        if self.witness.input_state_id != intermediate.post_temporal_state_id:
            raise ValueError("local rank solve input state is not the literal postTemporal state")
        if self.witness.debt_family_id != intermediate.debt_family_id:
            raise ValueError("local rank solve debt family is not the actual postTemporal debt")
        if self.witness.band_index != intermediate.band_index:
            raise ValueError("local rank solve band does not match the admitted debt evaluation")
        if self.witness.slow_point_id != intermediate.slow_point_id:
            raise ValueError("local rank solve slow point does not match the admitted debt evaluation")
        if rank.debt_id != intermediate.band_debt_id:
            raise ValueError("compact rank repair is not bound to the admitted band debt")

    @property
    def local_five_row_solve_admitted(self) -> bool:
        return True

    @property
    def linear_debt_cancellation_admitted(self) -> bool:
        return True

    @property
    def required_masses_preserved_admitted(self) -> bool:
        return True

    @property
    def post_rank_debt_equals_remainders_admitted(self) -> bool:
        return True

    @property
    def moving_support_preserved_admitted(self) -> bool:
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
    def actual_remainder_values_materialized(self) -> bool:
        return False

    @property
    def section9_iteration_input_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
