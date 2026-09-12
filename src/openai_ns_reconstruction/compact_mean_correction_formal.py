"""Fail-closed formal admission for the compact Section 8 mean correction.

The executable helper :mod:`mean_rank_update` is intentionally diagnostic: its
compact polynomial bump is *not* the noncomputable Mathlib bump used by the
pinned Lean development.  This module provides the parallel paper-exact handoff
instead.  It pins ``CorrectionState.rank_rows_on_patch``, whose formal statement
uses the actual correction-state debt and the physical five-row construction,
and only accepts an external ``lean-formal-export`` theorem application.

No Python bump, debt vector, background field, or corrected velocity is
materialized here.  In particular, this module does not identify the finite-head
signed mean defect with ``CorrectionState.debt``; that bridge needs its own
formal theorem/export and is deliberately left open.
"""
from __future__ import annotations

from dataclasses import dataclass


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/CorrectionState.lean"
PINNED_RANK_ROWS_THEOREM = "NavierStokes.CorrectionState.rank_rows_on_patch"
PINNED_FIVE_ROWS_DEPENDENCY = "NavierStokes.MeanRankUpdate.physical_five_rows"
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_RANK_ROWS_THEOREM,
    PINNED_FIVE_ROWS_DEPENDENCY,
)
PRODUCER_KIND = "lean-formal-export"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


@dataclass(frozen=True)
class CompactMeanCorrectionWitness:
    """Opaque metadata for one ``rank_rows_on_patch`` theorem application.

    ``debt_id`` names the exact ``CorrectionState.debt c u n x`` object used by
    Lean.  The two background-model identifiers name the angular and axial
    profiles whose agreement with the actual base is assumed on the support of
    the constructed bumps.  These are identities only; Python does not encode
    the corresponding functions.
    """

    rank_data_id: str
    context_id: str
    state_id: str
    debt_id: str
    slow_point_id: str
    plane_point_id: str
    angular_background_id: str
    axial_background_id: str
    application_id: str
    producer_kind: str
    provenance: str
    band_index: int

    lambda_positive_certified: bool
    coefficient_nonzero_certified: bool
    inner_positive_certified: bool
    inner_lt_outer_certified: bool
    length_positive_certified: bool
    velocity_nonzero_certified: bool
    angular_background_agreement_on_support_certified: bool
    axial_background_agreement_on_support_certified: bool
    physical_five_rows_dependency_certified: bool
    theorem_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_RANK_ROWS_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "rank_data_id",
            "context_id",
            "state_id",
            "debt_id",
            "slow_point_id",
            "plane_point_id",
            "angular_background_id",
            "axial_background_id",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if not isinstance(self.band_index, int) or isinstance(self.band_index, bool):
            raise TypeError("band_index must be an integer")
        if self.band_index < 0:
            raise ValueError("band_index must be nonnegative")

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must match CorrectionState.lean exactly")
        if self.theorem_symbol != PINNED_RANK_ROWS_THEOREM:
            raise ValueError("theorem_symbol must be CorrectionState.rank_rows_on_patch")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned five-row chain exactly")

        for name in (
            "lambda_positive_certified",
            "coefficient_nonzero_certified",
            "inner_positive_certified",
            "inner_lt_outer_certified",
            "length_positive_certified",
            "velocity_nonzero_certified",
            "angular_background_agreement_on_support_certified",
            "axial_background_agreement_on_support_certified",
            "physical_five_rows_dependency_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class CompactMeanCorrectionAdmission:
    """Validated formal handoff for one compact five-row mean correction."""

    witness: CompactMeanCorrectionWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, CompactMeanCorrectionWitness):
            raise TypeError("witness must be a CompactMeanCorrectionWitness")

    @property
    def compact_five_row_theorem_admitted(self) -> bool:
        return True

    @property
    def actual_state_debt_identity_bound(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def theorem_application_machine_replayed(self) -> bool:
        return False

    @property
    def actual_debt_values_materialized(self) -> bool:
        return False

    @property
    def actual_bump_values_materialized(self) -> bool:
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
