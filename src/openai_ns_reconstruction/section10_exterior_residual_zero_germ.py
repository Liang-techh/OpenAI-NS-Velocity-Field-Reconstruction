"""Fail-closed theorem-shaped certificate for exterior residual zero germs.

The pinned formal proof of ``LocalPaper.selected_exterior_residual_zero`` is
stronger than its value-level conclusion at one point.  For every point in the
strict paper exterior, ``LocalPaperDomain.exterior_conditions_eventually``
keeps all three strict exterior hypotheses true on a neighborhood.  Applying
``selected_exterior_residual_zero`` pointwise on that neighborhood gives local
equality of the literal residual with zero.  The pinned theorem
``SolenoidalDiagonal.iteratedFDeriv_eventuallyEq`` then transports that zero
germ to every actual Frechet derivative order.

This module records that exact proof chain for the *same* schedule and local
velocity/pressure already admitted by the paper-localization spine.  It does
not replay Lean, materialize the fields, construct a force, or turn a formal
zero-germ consequence into runtime support verification.
"""
from __future__ import annotations

from dataclasses import dataclass

from openai_ns_reconstruction.section10_paper_localization_spine import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
    Section10PaperLocalizationSpineAdmission,
)


PINNED_LOCAL_PAPER_FILE = "NavierStokes/LocalPaperTheorem.lean"
PINNED_LOCAL_DOMAIN_FILE = "NavierStokes/LocalPaperDomain.lean"
PINNED_DERIVATIVE_FILE = "NavierStokes/SolenoidalDiagonal.lean"
PINNED_EXTERIOR_POINT_THEOREM = "NavierStokes.LocalPaper.selected_exterior_residual_zero"
PINNED_EXTERIOR_NEIGHBORHOOD_THEOREM = (
    "NavierStokes.LocalPaperDomain.exterior_conditions_eventually"
)
PINNED_DERIVATIVE_GERM_THEOREM = "NavierStokes.SolenoidalDiagonal.iteratedFDeriv_eventuallyEq"
PINNED_ZERO_DERIVATIVE_IDENTITY = "iteratedFDeriv_fun_zero"
PINNED_EVIDENCE_KIND = "lean-proof-chain-export"
UNIVERSAL_DERIVATIVE_QUANTIFIER = "forall-m-in-N"
ZERO_GERM_KIND = "nhds-eventual-equality-to-zero"
STRICT_EXTERIOR_SCOPE = "LocalPaperDomain.strict-exterior-neighborhood"

REQUIRED_PROOF_SYMBOLS = (
    PINNED_EXTERIOR_POINT_THEOREM,
    PINNED_EXTERIOR_NEIGHBORHOOD_THEOREM,
    PINNED_DERIVATIVE_GERM_THEOREM,
    PINNED_ZERO_DERIVATIVE_IDENTITY,
)
REQUIRED_FORMAL_FILES = (
    PINNED_LOCAL_PAPER_FILE,
    PINNED_LOCAL_DOMAIN_FILE,
    PINNED_DERIVATIVE_FILE,
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
class Section10ExteriorResidualZeroGermWitness:
    """Exact proof-chain metadata for one paper-local exterior residual germ."""

    application_id: str
    schedule_id: str
    local_velocity_id: str
    local_pressure_id: str
    residual_id: str
    provenance: str

    evidence_kind: str = PINNED_EVIDENCE_KIND
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_files: tuple[str, ...] = REQUIRED_FORMAL_FILES
    exterior_point_theorem_symbol: str = PINNED_EXTERIOR_POINT_THEOREM
    exterior_neighborhood_theorem_symbol: str = PINNED_EXTERIOR_NEIGHBORHOOD_THEOREM
    derivative_germ_theorem_symbol: str = PINNED_DERIVATIVE_GERM_THEOREM
    zero_derivative_identity_symbol: str = PINNED_ZERO_DERIVATIVE_IDENTITY
    derivative_order_quantifier: str = UNIVERSAL_DERIVATIVE_QUANTIFIER
    zero_germ_kind: str = ZERO_GERM_KIND
    exterior_scope: str = STRICT_EXTERIOR_SCOPE
    dependency_symbols: tuple[str, ...] = REQUIRED_PROOF_SYMBOLS

    max_derivative_order: int | None = None
    sampled_or_fitted_evidence: bool = False
    manufactured_residual: bool = False
    force_defined_as_residual: bool = False

    same_schedule_as_paper_spine_certified: bool = False
    same_local_velocity_as_paper_spine_certified: bool = False
    same_local_pressure_as_paper_spine_certified: bool = False
    literal_navier_stokes_residual_certified: bool = False
    strict_exterior_conditions_certified: bool = False
    exterior_conditions_persist_on_neighborhood_certified: bool = False
    pointwise_residual_zero_on_persisting_neighborhood_certified: bool = False
    residual_zero_germ_certified: bool = False
    all_order_exterior_residual_jets_zero_certified: bool = False

    def __post_init__(self) -> None:
        for name in (
            "application_id",
            "schedule_id",
            "local_velocity_id",
            "local_pressure_id",
            "residual_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        expected = {
            "evidence_kind": PINNED_EVIDENCE_KIND,
            "formal_repository": PINNED_FORMAL_REPOSITORY,
            "formal_commit": PINNED_FORMAL_COMMIT,
            "formal_files": REQUIRED_FORMAL_FILES,
            "exterior_point_theorem_symbol": PINNED_EXTERIOR_POINT_THEOREM,
            "exterior_neighborhood_theorem_symbol": PINNED_EXTERIOR_NEIGHBORHOOD_THEOREM,
            "derivative_germ_theorem_symbol": PINNED_DERIVATIVE_GERM_THEOREM,
            "zero_derivative_identity_symbol": PINNED_ZERO_DERIVATIVE_IDENTITY,
            "derivative_order_quantifier": UNIVERSAL_DERIVATIVE_QUANTIFIER,
            "zero_germ_kind": ZERO_GERM_KIND,
            "exterior_scope": STRICT_EXTERIOR_SCOPE,
            "dependency_symbols": REQUIRED_PROOF_SYMBOLS,
        }
        for name, value in expected.items():
            if getattr(self, name) != value:
                raise ValueError(f"{name} must match the pinned exterior residual proof chain exactly")

        if self.max_derivative_order is not None:
            raise ValueError("finite derivative frontiers cannot certify all-order exterior zero jets")
        if self.sampled_or_fitted_evidence is not False:
            raise ValueError("sampled/fitted evidence cannot certify a residual zero germ")
        if self.manufactured_residual is not False:
            raise ValueError("manufactured residual evidence is forbidden")
        if self.force_defined_as_residual is not False:
            raise ValueError("f=R is not independent force verification")

        for name in (
            "same_schedule_as_paper_spine_certified",
            "same_local_velocity_as_paper_spine_certified",
            "same_local_pressure_as_paper_spine_certified",
            "literal_navier_stokes_residual_certified",
            "strict_exterior_conditions_certified",
            "exterior_conditions_persist_on_neighborhood_certified",
            "pointwise_residual_zero_on_persisting_neighborhood_certified",
            "residual_zero_germ_certified",
            "all_order_exterior_residual_jets_zero_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10ExteriorResidualZeroGermAdmission:
    """Formal all-order exterior zero-jet binding; runtime flags stay false."""

    spine: Section10PaperLocalizationSpineAdmission
    witness: Section10ExteriorResidualZeroGermWitness

    def __post_init__(self) -> None:
        if not isinstance(self.spine, Section10PaperLocalizationSpineAdmission):
            raise TypeError("spine must be a Section10PaperLocalizationSpineAdmission")
        if not isinstance(self.witness, Section10ExteriorResidualZeroGermWitness):
            raise TypeError("witness must be a Section10ExteriorResidualZeroGermWitness")

        parent = self.spine.witness
        for child_name, parent_name in (
            ("schedule_id", "schedule_id"),
            ("local_velocity_id", "local_velocity_id"),
            ("local_pressure_id", "local_pressure_id"),
        ):
            if getattr(self.witness, child_name) != getattr(parent, parent_name):
                raise ValueError(
                    f"{child_name} must be identical to the admitted paper-localization spine"
                )

    @property
    def formal_exterior_residual_zero_germ_admitted(self) -> bool:
        return True

    @property
    def formal_all_order_exterior_residual_jets_zero_admitted(self) -> bool:
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
    def support_exterior_zero_jets_runtime_verified(self) -> bool:
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
    def paper_exact_velocity_available(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False


def admit_exterior_residual_zero_germ(
    spine: Section10PaperLocalizationSpineAdmission,
    export: Section10ExteriorResidualZeroGermWitness,
) -> Section10ExteriorResidualZeroGermAdmission:
    """Bind the exact all-order exterior zero-germ proof chain to one spine."""

    return Section10ExteriorResidualZeroGermAdmission(spine=spine, witness=export)
