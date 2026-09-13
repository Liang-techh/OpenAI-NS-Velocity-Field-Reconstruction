"""Fail-closed binding from closed actual stage estimates to the selected schedule.

``LocalResidualFlatness.selectedEstimates`` is not an independent estimate
oracle.  At the pinned formal revision it is definitionally the closed
specialization of ``ActualCandidateAssembly.estimates`` at the repository's
``selectedBudget``/``selectedThreshold`` choice, and that candidate estimate is
assembled through ``GluedStageEstimates.actualStageEstimates`` from the actual
stage-estimate construction.

The already-landed selected-schedule admission names ``selectedEstimates`` and
its selected potential/direct/pressure families, but by itself an exported
metadata record could still be cross-wired to unrelated opaque identities.
This module closes exactly that identity seam.  It does not evaluate a stage,
choose the existential schedule, replay Lean, or materialize the paper field.
"""
from __future__ import annotations

from dataclasses import dataclass

from .section9_selected_schedule_flatness import (
    Section9SelectedScheduleFlatnessAdmission,
)


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILES = (
    "NavierStokes/ActualCandidateAssembly.lean",
    "NavierStokes/LocalResidualFlatness.lean",
    "NavierStokes/ActualStageEstimates.lean",
)
PINNED_SELECTED_ESTIMATES = "NavierStokes.LocalResidualFlatness.selectedEstimates"
PINNED_CANDIDATE_ESTIMATES = "NavierStokes.ActualCandidateAssembly.estimates"
PINNED_GLUED_STAGE_ESTIMATES = "NavierStokes.GluedStageEstimates.actualStageEstimates"
PINNED_STAGE_ESTIMATES = (
    "NavierStokes.ActualStageEstimates.stageEstimates_of_representations"
)
PINNED_SELECTED_BUDGET = "NavierStokes.ActualCandidateConstruction.selectedBudget"
PINNED_SELECTED_THRESHOLD = "NavierStokes.ActualCandidateConstruction.selectedThreshold"
PINNED_SELECTED_THRESHOLD_GEOMETRY = (
    "NavierStokes.ActualCandidateConstruction.selectedThreshold_geometry"
)
PINNED_SELECTED_QBIG = "NavierStokes.ActualCandidateConstruction.selectedQbig"
PINNED_SELECTED_POTENTIAL_STAGES = (
    "NavierStokes.ActualCandidateAssembly.selectedPotentialStages"
)
PINNED_SELECTED_DIRECT_STAGES = (
    "NavierStokes.ActualCandidateAssembly.selectedDirectStages"
)
PINNED_SELECTED_PRESSURE_STAGES = (
    "NavierStokes.ActualCandidateAssembly.selectedPressureStages"
)
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_SELECTED_ESTIMATES,
    PINNED_CANDIDATE_ESTIMATES,
    PINNED_GLUED_STAGE_ESTIMATES,
    PINNED_STAGE_ESTIMATES,
    PINNED_SELECTED_BUDGET,
    PINNED_SELECTED_THRESHOLD,
    PINNED_SELECTED_THRESHOLD_GEOMETRY,
    PINNED_SELECTED_QBIG,
    PINNED_SELECTED_POTENTIAL_STAGES,
    PINNED_SELECTED_DIRECT_STAGES,
    PINNED_SELECTED_PRESSURE_STAGES,
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
class Section9SelectedEstimatesBindingWitness:
    """Opaque identities for the exact closed specialization used by Section 9."""

    selected_budget_id: str
    selected_threshold_id: str
    selected_threshold_geometry_id: str
    selected_qbig_id: str
    candidate_estimates_id: str
    selected_estimates_id: str
    potential_stages_id: str
    direct_stages_id: str
    pressure_stages_id: str
    candidate_estimates_application_id: str
    binding_application_id: str
    producer_kind: str
    provenance: str

    selected_parameters_identified_certified: bool
    selected_threshold_geometry_certified: bool
    selected_qbig_specialization_certified: bool
    selected_stage_families_specialization_certified: bool
    candidate_estimates_actual_constructor_certified: bool
    selected_estimates_eq_candidate_certified: bool
    binding_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_files: tuple[str, ...] = PINNED_FORMAL_FILES
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "selected_budget_id",
            "selected_threshold_id",
            "selected_threshold_geometry_id",
            "selected_qbig_id",
            "candidate_estimates_id",
            "selected_estimates_id",
            "potential_stages_id",
            "direct_stages_id",
            "pressure_stages_id",
            "candidate_estimates_application_id",
            "binding_application_id",
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
        if self.formal_files != PINNED_FORMAL_FILES:
            raise ValueError("formal_files must match the pinned closed-estimate chain exactly")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned closed-estimate chain exactly")

        for name in (
            "selected_parameters_identified_certified",
            "selected_threshold_geometry_certified",
            "selected_qbig_specialization_certified",
            "selected_stage_families_specialization_certified",
            "candidate_estimates_actual_constructor_certified",
            "selected_estimates_eq_candidate_certified",
            "binding_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section9SelectedEstimatesBindingAdmission:
    """Bind ``selected_schedule`` to the exact closed actual estimate record."""

    selected_schedule: Section9SelectedScheduleFlatnessAdmission
    witness: Section9SelectedEstimatesBindingWitness

    def __post_init__(self) -> None:
        if not isinstance(
            self.selected_schedule, Section9SelectedScheduleFlatnessAdmission
        ):
            raise TypeError(
                "selected_schedule must be a Section9SelectedScheduleFlatnessAdmission"
            )
        if not isinstance(self.witness, Section9SelectedEstimatesBindingWitness):
            raise TypeError("witness must be a Section9SelectedEstimatesBindingWitness")

        schedule = self.selected_schedule.witness
        pairs = (
            ("selected_estimates_id", self.witness.selected_estimates_id, schedule.selected_estimates_id),
            ("selected_qbig_id", self.witness.selected_qbig_id, schedule.selected_qbig_id),
            ("potential_stages_id", self.witness.potential_stages_id, schedule.potential_stages_id),
            ("direct_stages_id", self.witness.direct_stages_id, schedule.direct_stages_id),
            ("pressure_stages_id", self.witness.pressure_stages_id, schedule.pressure_stages_id),
        )
        for name, closed_id, schedule_id in pairs:
            if closed_id != schedule_id:
                raise ValueError(
                    f"{name} must be the identical closed object used by selected_schedule"
                )

        if schedule.formal_repository != self.witness.formal_repository:
            raise ValueError("formal_repository mismatch across the composed admissions")
        if schedule.formal_commit != self.witness.formal_commit:
            raise ValueError("formal_commit mismatch across the composed admissions")

    @property
    def selected_schedule_bound_to_closed_actual_estimates(self) -> bool:
        return True

    @property
    def selected_stage_families_bound_to_closed_actual_choice(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def closed_estimate_binding_machine_replayed(self) -> bool:
        return False

    @property
    def selected_estimates_values_materialized(self) -> bool:
        return False

    @property
    def selected_stage_fields_materialized(self) -> bool:
        return False

    @property
    def selected_schedule_values_materialized(self) -> bool:
        return False

    @property
    def section9_iteration_machine_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
