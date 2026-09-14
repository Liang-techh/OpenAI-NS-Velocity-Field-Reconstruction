"""Fail-closed handoff for the closed actual all-order residual-flatness theorem.

At the pinned OpenAI formal revision, ``LocalResidualFlatness.selected_schedule``
is a *closed* theorem for the repository's selected actual stage families.  It
chooses one schedule once and returns the same schedule together with
``ThreeCutBounds`` and ``AllResidualJetRates``.  The latter quantifies over every
natural derivative order ``m`` and every nonnegative real decay order ``r`` for
the literal ``MixedDiagonalResidual.residual`` at the open-past endpoint
``(t,x) -> (1,0)``.

This is materially stronger than a finite derivative ladder or a collection of
sampled rates: the existential schedule is outside both order quantifiers.  It
also removes the raw ``StageEstimates`` hypotheses at this interface because
``selected_schedule`` internally uses ``selectedEstimates`` from the actual
fixed candidate assembly.

Scope is intentionally narrow.  The theorem is local at the singular spacetime
point and concerns the unlocalized mixed residual.  It does *not* by itself
materialize Python fields, prove uniform bounds on the entire official
``3/4 <= t < 1`` plateau, prove support-exterior zero neighborhoods, or construct
the Section 10 force.  Those remain separate closure obligations.
"""
from __future__ import annotations

from dataclasses import dataclass


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/LocalResidualFlatness.lean"
PINNED_THEOREM = "NavierStokes.LocalResidualFlatness.selected_schedule"
PINNED_RATE_PREDICATE = "NavierStokes.LocalResidualFlatness.AllResidualJetRates"
PINNED_SELECTED_ESTIMATES = "NavierStokes.LocalResidualFlatness.selectedEstimates"
PINNED_POTENTIAL_STAGES = "NavierStokes.ActualCandidateAssembly.selectedPotentialStages"
PINNED_DIRECT_STAGES = "NavierStokes.ActualCandidateAssembly.selectedDirectStages"
PINNED_PRESSURE_STAGES = "NavierStokes.ActualCandidateAssembly.selectedPressureStages"
PINNED_RESIDUAL = "NavierStokes.MixedDiagonalResidual.residual"
PINNED_PHYSICAL_SCALE = "NavierStokes.PhysicalWaveSum.physicalQ"
PINNED_ENDPOINT_FILTER = "nhdsWithin-(1,0)-openPast-1"
PRODUCER_KIND = "lean-formal-export"
UNIVERSAL_DERIVATIVE_QUANTIFIER = "forall-m-in-N"
UNIVERSAL_DECAY_QUANTIFIER = "forall-r-in-R-with-0<=r"

REQUIRED_FORMAL_SYMBOLS = (
    PINNED_THEOREM,
    PINNED_RATE_PREDICATE,
    PINNED_SELECTED_ESTIMATES,
    PINNED_POTENTIAL_STAGES,
    PINNED_DIRECT_STAGES,
    PINNED_PRESSURE_STAGES,
    PINNED_RESIDUAL,
    PINNED_PHYSICAL_SCALE,
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
class Section9ActualLocalResidualFlatnessWitness:
    """Opaque export metadata for one application of the closed actual theorem."""

    schedule_id: str
    application_id: str
    provenance: str

    producer_kind: str = PRODUCER_KIND
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_THEOREM
    rate_predicate_symbol: str = PINNED_RATE_PREDICATE
    selected_estimates_id: str = PINNED_SELECTED_ESTIMATES
    potential_stage_family_id: str = PINNED_POTENTIAL_STAGES
    direct_stage_family_id: str = PINNED_DIRECT_STAGES
    pressure_stage_family_id: str = PINNED_PRESSURE_STAGES
    residual_id: str = PINNED_RESIDUAL
    physical_scale_id: str = PINNED_PHYSICAL_SCALE
    endpoint_filter_id: str = PINNED_ENDPOINT_FILTER
    derivative_order_quantifier: str = UNIVERSAL_DERIVATIVE_QUANTIFIER
    decay_order_quantifier: str = UNIVERSAL_DECAY_QUANTIFIER
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    max_derivative_order: int | None = None
    max_decay_order: float | None = None
    sampled_or_fitted_evidence: bool = False
    manufactured_residual: bool = False

    closed_actual_theorem_application_certified: bool = False
    one_schedule_outside_order_quantifiers_certified: bool = False
    selected_schedule_certified: bool = False
    three_cut_bounds_certified: bool = False
    all_residual_jet_rates_certified: bool = False
    actual_unlocalized_mixed_residual_certified: bool = False
    open_past_endpoint_filter_certified: bool = False

    def __post_init__(self) -> None:
        for name in ("schedule_id", "application_id", "provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        expected = {
            "producer_kind": PRODUCER_KIND,
            "formal_repository": PINNED_FORMAL_REPOSITORY,
            "formal_commit": PINNED_FORMAL_COMMIT,
            "formal_file": PINNED_FORMAL_FILE,
            "theorem_symbol": PINNED_THEOREM,
            "rate_predicate_symbol": PINNED_RATE_PREDICATE,
            "selected_estimates_id": PINNED_SELECTED_ESTIMATES,
            "potential_stage_family_id": PINNED_POTENTIAL_STAGES,
            "direct_stage_family_id": PINNED_DIRECT_STAGES,
            "pressure_stage_family_id": PINNED_PRESSURE_STAGES,
            "residual_id": PINNED_RESIDUAL,
            "physical_scale_id": PINNED_PHYSICAL_SCALE,
            "endpoint_filter_id": PINNED_ENDPOINT_FILTER,
            "derivative_order_quantifier": UNIVERSAL_DERIVATIVE_QUANTIFIER,
            "decay_order_quantifier": UNIVERSAL_DECAY_QUANTIFIER,
            "dependency_symbols": REQUIRED_FORMAL_SYMBOLS,
        }
        for name, value in expected.items():
            if getattr(self, name) != value:
                raise ValueError(f"{name} must match the pinned actual residual-flatness theorem exactly")

        if self.max_derivative_order is not None:
            raise ValueError("finite derivative frontiers are not all-order residual evidence")
        if self.max_decay_order is not None:
            raise ValueError("finite decay-order frontiers are not all-order residual evidence")
        if self.sampled_or_fitted_evidence is not False:
            raise ValueError("sampled/fitted evidence cannot certify the closed formal theorem")
        if self.manufactured_residual is not False:
            raise ValueError("manufactured residual evidence is forbidden")

        for name in (
            "closed_actual_theorem_application_certified",
            "one_schedule_outside_order_quantifiers_certified",
            "selected_schedule_certified",
            "three_cut_bounds_certified",
            "all_residual_jet_rates_certified",
            "actual_unlocalized_mixed_residual_certified",
            "open_past_endpoint_filter_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section9ActualLocalResidualFlatnessAdmission:
    """Validated theorem-facing admission; never a runtime field artifact."""

    witness: Section9ActualLocalResidualFlatnessWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section9ActualLocalResidualFlatnessWitness):
            raise TypeError("witness must be a Section9ActualLocalResidualFlatnessWitness")

    @property
    def actual_selected_schedule_all_order_flatness_admitted(self) -> bool:
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
    def official_late_plateau_uniform_bounds_verified(self) -> bool:
        return False

    @property
    def support_exterior_all_order_zero_verified(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def forcing_artifact_ready(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False


def admit_actual_local_residual_flatness(
    export: Section9ActualLocalResidualFlatnessWitness,
) -> Section9ActualLocalResidualFlatnessAdmission:
    """Admit only the exact closed formal export after dataclass validation."""

    if not isinstance(export, Section9ActualLocalResidualFlatnessWitness):
        raise TypeError("export must be a Section9ActualLocalResidualFlatnessWitness")
    return Section9ActualLocalResidualFlatnessAdmission(export)
