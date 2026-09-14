"""Fail-closed handoff for the paper-specific local-to-compact theorem spine.

The pinned formal repository contains the closed theorem
``PaperLocalization.local_theorem_with_compact_candidate``.  Its proof chooses
one actual diagonal schedule using ``LocalResidualFlatness.selected_schedule``
and uses that *same* schedule for both the local paper fields and the compact
whole-space candidate.  The theorem simultaneously returns:

* ``LocalPaper.Properties`` for the local fields, including all-order residual
  flatness, genuine away-from-origin endpoint extensions, divergence-free flow,
  and the literal exterior residual-zero statement;
* ``NavierStokesR3.ProblemStatement.CandidateProperties 1`` for the localized
  whole-space velocity/pressure/force, including smooth compact positive-time
  forcing, exact forced Navier--Stokes, divergence-free flow, one uniform
  kinetic-energy bound on ``0 <= t < 1``, and speed blow-up at ``t=1``;
* nonexistence of a global smooth finite-energy solution for the same force and
  zero datum; and
* exact local/compact field agreement on one open neighborhood of the origin
  for every official late time ``3/4 <= t < 1``.

This module records only that exact formal theorem handoff.  It does not replay
Lean and it does not materialize the theorem's existential fields in Python.
In particular, formal existence and same-field agreement are not runtime
residual/force artifacts or an independent numerical reconstruction.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/PaperLocalization.lean"
PINNED_THEOREM = "NavierStokes.PaperLocalization.local_theorem_with_compact_candidate"
PINNED_SELECTED_SCHEDULE = "NavierStokes.LocalResidualFlatness.selected_schedule"
PINNED_LOCAL_PROPERTIES_OF_SCHEDULE = "NavierStokes.LocalPaper.properties_of_schedule"
PINNED_SELECTED_COMPACT_CANDIDATE = (
    "NavierStokes.LocalScheduleWitness.selected_compact_candidate"
)
PINNED_LOCAL_PROPERTIES = "NavierStokes.LocalPaper.Properties"
PINNED_CANDIDATE_PROPERTIES = "NavierStokesR3.ProblemStatement.CandidateProperties"
PINNED_PLATEAU = "NavierStokes.SpatialLocalization.plateau"
PINNED_EVIDENCE_KIND = "lean-formal-export"

UNIT_VISCOSITY = Fraction(1, 1)
OFFICIAL_LATE_START = Fraction(3, 4)
OFFICIAL_ENDPOINT = Fraction(1, 1)
UNIVERSAL_DERIVATIVE_QUANTIFIER = "forall-m-in-N"
UNIVERSAL_DECAY_QUANTIFIER = "forall-r-in-R-with-0<=r"

REQUIRED_FORMAL_SYMBOLS = (
    PINNED_THEOREM,
    PINNED_SELECTED_SCHEDULE,
    PINNED_LOCAL_PROPERTIES_OF_SCHEDULE,
    PINNED_SELECTED_COMPACT_CANDIDATE,
    PINNED_LOCAL_PROPERTIES,
    PINNED_CANDIDATE_PROPERTIES,
    PINNED_PLATEAU,
)


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


@dataclass(frozen=True)
class Section10PaperLocalizationSpineWitness:
    """Opaque identifiers for one exact application of the closed paper theorem."""

    application_id: str
    schedule_id: str
    local_potential_id: str
    local_direct_id: str
    local_velocity_id: str
    local_pressure_id: str
    compact_velocity_id: str
    compact_pressure_id: str
    force_id: str
    compact_support_id: str
    plateau_open_set_id: str
    provenance: str

    evidence_kind: str = PINNED_EVIDENCE_KIND
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_THEOREM
    selected_schedule_symbol: str = PINNED_SELECTED_SCHEDULE
    local_properties_of_schedule_symbol: str = PINNED_LOCAL_PROPERTIES_OF_SCHEDULE
    selected_compact_candidate_symbol: str = PINNED_SELECTED_COMPACT_CANDIDATE
    local_properties_symbol: str = PINNED_LOCAL_PROPERTIES
    candidate_properties_symbol: str = PINNED_CANDIDATE_PROPERTIES
    plateau_set_id: str = PINNED_PLATEAU
    viscosity: Fraction = UNIT_VISCOSITY
    late_start: Fraction = OFFICIAL_LATE_START
    endpoint: Fraction = OFFICIAL_ENDPOINT
    derivative_order_quantifier: str = UNIVERSAL_DERIVATIVE_QUANTIFIER
    decay_order_quantifier: str = UNIVERSAL_DECAY_QUANTIFIER
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    max_derivative_order: int | None = None
    sampled_or_fitted_evidence: bool = False
    manufactured_residual: bool = False
    generic_localization_substitute: bool = False

    closed_theorem_application_certified: bool = False
    same_selected_schedule_certified: bool = False
    local_properties_certified: bool = False
    local_away_extensions_certified: bool = False
    local_all_order_residual_flatness_certified: bool = False
    local_exterior_residual_zero_certified: bool = False
    compact_candidate_properties_certified: bool = False
    compact_force_smooth_certified: bool = False
    compact_force_positive_time_support_certified: bool = False
    compact_spatial_support_certified: bool = False
    compact_divergence_free_certified: bool = False
    compact_exact_navier_stokes_certified: bool = False
    compact_uniform_finite_energy_certified: bool = False
    compact_speed_unbounded_certified: bool = False
    no_global_finite_energy_solution_certified: bool = False
    initial_rest_certified: bool = False
    plateau_open_neighborhood_contains_origin_certified: bool = False
    official_late_plateau_same_fields_certified: bool = False

    def __post_init__(self) -> None:
        for name in (
            "application_id",
            "schedule_id",
            "local_potential_id",
            "local_direct_id",
            "local_velocity_id",
            "local_pressure_id",
            "compact_velocity_id",
            "compact_pressure_id",
            "force_id",
            "compact_support_id",
            "plateau_open_set_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        expected = {
            "evidence_kind": PINNED_EVIDENCE_KIND,
            "formal_repository": PINNED_FORMAL_REPOSITORY,
            "formal_commit": PINNED_FORMAL_COMMIT,
            "formal_file": PINNED_FORMAL_FILE,
            "theorem_symbol": PINNED_THEOREM,
            "selected_schedule_symbol": PINNED_SELECTED_SCHEDULE,
            "local_properties_of_schedule_symbol": PINNED_LOCAL_PROPERTIES_OF_SCHEDULE,
            "selected_compact_candidate_symbol": PINNED_SELECTED_COMPACT_CANDIDATE,
            "local_properties_symbol": PINNED_LOCAL_PROPERTIES,
            "candidate_properties_symbol": PINNED_CANDIDATE_PROPERTIES,
            "plateau_set_id": PINNED_PLATEAU,
            "derivative_order_quantifier": UNIVERSAL_DERIVATIVE_QUANTIFIER,
            "decay_order_quantifier": UNIVERSAL_DECAY_QUANTIFIER,
            "dependency_symbols": REQUIRED_FORMAL_SYMBOLS,
        }
        for name, value in expected.items():
            if getattr(self, name) != value:
                raise ValueError(f"{name} must match the pinned paper-localization theorem exactly")

        for name, expected_fraction in (
            ("viscosity", UNIT_VISCOSITY),
            ("late_start", OFFICIAL_LATE_START),
            ("endpoint", OFFICIAL_ENDPOINT),
        ):
            value = getattr(self, name)
            if type(value) is not Fraction or value != expected_fraction:
                raise ValueError(f"{name} must be the exact pinned rational {expected_fraction}")

        if self.max_derivative_order is not None:
            raise ValueError("finite derivative frontiers cannot certify the all-order local theorem")
        if self.sampled_or_fitted_evidence is not False:
            raise ValueError("sampled/fitted evidence cannot certify the closed formal theorem")
        if self.manufactured_residual is not False:
            raise ValueError("manufactured residual evidence is forbidden")
        if self.generic_localization_substitute is not False:
            raise ValueError("generic localization cannot substitute for PaperLocalization")

        for name in (
            "closed_theorem_application_certified",
            "same_selected_schedule_certified",
            "local_properties_certified",
            "local_away_extensions_certified",
            "local_all_order_residual_flatness_certified",
            "local_exterior_residual_zero_certified",
            "compact_candidate_properties_certified",
            "compact_force_smooth_certified",
            "compact_force_positive_time_support_certified",
            "compact_spatial_support_certified",
            "compact_divergence_free_certified",
            "compact_exact_navier_stokes_certified",
            "compact_uniform_finite_energy_certified",
            "compact_speed_unbounded_certified",
            "no_global_finite_energy_solution_certified",
            "initial_rest_certified",
            "plateau_open_neighborhood_contains_origin_certified",
            "official_late_plateau_same_fields_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10PaperLocalizationSpineAdmission:
    """Formal end-to-end theorem binding; runtime reconstruction stays false."""

    witness: Section10PaperLocalizationSpineWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section10PaperLocalizationSpineWitness):
            raise TypeError("witness must be a Section10PaperLocalizationSpineWitness")

    @property
    def formal_paper_localization_spine_admitted(self) -> bool:
        return True

    @property
    def official_late_plateau_field_agreement_formally_bound(self) -> bool:
        return True

    @property
    def final_candidate_properties_formally_bound(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_fields_materialized(self) -> bool:
        return False

    @property
    def actual_section9_sequence_verified(self) -> bool:
        return False

    @property
    def section9_all_order_endpoint_limits_verified(self) -> bool:
        return False

    @property
    def residual_artifact_ready(self) -> bool:
        return False

    @property
    def forcing_artifact_ready(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def compact_support_runtime_verified(self) -> bool:
        return False

    @property
    def divergence_free_closure_verified(self) -> bool:
        return False

    @property
    def finite_energy_closure_verified(self) -> bool:
        return False

    @property
    def blow_up_closure_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False


def admit_paper_localization_spine(
    export: Section10PaperLocalizationSpineWitness,
) -> Section10PaperLocalizationSpineAdmission:
    """Admit only the exact closed PaperLocalization theorem export."""

    if not isinstance(export, Section10PaperLocalizationSpineWitness):
        raise TypeError("export must be a Section10PaperLocalizationSpineWitness")
    return Section10PaperLocalizationSpineAdmission(export)
