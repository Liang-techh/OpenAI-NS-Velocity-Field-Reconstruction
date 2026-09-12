"""Bind the exact actual signed mean-cross defect to the canonical Section 6/7 scope.

Section 8 cannot treat the finite low-band mean error as an arbitrary numerical
residual.  At the pinned formal revision,
``ActualSignedMeanBinding.requested_cross_defect`` proves the structural identity

``meanBar(actualCross) - requestedStress
   = -missingWeight(Prepared.N, physicalScale) * requestedStress``.

Unlike the tail theorem, this identity itself has no ``Prepared.N + 1 <= n``
hypothesis: it exposes the cutoff deficit on every admitted band.  The already
landed tail admission separately certifies exact cancellation once the missing
weight vanishes on the tail.

This module is deliberately only a fail-closed handoff for an external
``lean-formal-export`` theorem application.  Lean context/state/point values and
the formal expressions for the mean cross, requested stress, physical scale,
and missing weight remain opaque stable identifiers.  Python neither samples
nor evaluates them and does not synthesize a compact mean correction.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .actual_signed_canonical_scope import (
    LargeBandCanonicalSignedCommonCurlScopeAdmission,
)
from .actual_signed_common_curl import PINNED_LEAN_COMMIT, PINNED_LEAN_REPOSITORY


PINNED_MEAN_FILE = "NavierStokes/ActualSignedMeanBinding.lean"
PINNED_DEFECT_THEOREM = "NavierStokes.ActualSignedMeanBinding.requested_cross_defect"
PINNED_FACTOR_THEOREM = "NavierStokes.ActualSignedMeanBinding.requested_cross_factor"
PINNED_PARTITION_DEFECT_THEOREM = (
    "NavierStokes.ActualPrimaryCovariance.partitionFactor_eq_one_sub_missing"
)
REQUIRED_DEFECT_SYMBOLS = (
    PINNED_DEFECT_THEOREM,
    PINNED_FACTOR_THEOREM,
    PINNED_PARTITION_DEFECT_THEOREM,
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
class ActualSignedMeanDefectWitness:
    """Opaque formal-export metadata for one ``requested_cross_defect`` use."""

    B: int
    N0: int
    prepared_N: int
    band: int
    component: int
    cycle_context_repr: str
    cycle_state_repr: str
    point_repr: str
    mean_cross_repr: str
    requested_stress_repr: str
    physical_scale_repr: str
    missing_weight_repr: str
    application_id: str
    producer_kind: str
    provenance: str
    canonical_active_family_identity_certified: bool
    actual_strip_membership_certified: bool
    actual_cross_identity_certified: bool
    requested_stress_identity_certified: bool
    physical_scale_identity_certified: bool
    prepared_N_identity_certified: bool
    theorem_application_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_MEAN_FILE
    theorem_symbol: str = PINNED_DEFECT_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_DEFECT_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "prepared_N", _natural(self.prepared_N, "prepared_N"))
        object.__setattr__(self, "band", _natural(self.band, "band"))
        component = _natural(self.component, "component")
        if component not in (0, 1):
            raise ValueError("component must be 0 (theta) or 1 (axial)")
        object.__setattr__(self, "component", component)

        for name in (
            "cycle_context_repr",
            "cycle_state_repr",
            "point_repr",
            "mean_cross_repr",
            "requested_stress_repr",
            "physical_scale_repr",
            "missing_weight_repr",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.lean_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("lean_repository must match the pinned formal source exactly")
        if self.lean_commit != PINNED_LEAN_COMMIT:
            raise ValueError("lean_commit must match the pinned formal source exactly")
        if self.lean_file != PINNED_MEAN_FILE:
            raise ValueError("lean_file must match the pinned actual signed mean source exactly")
        if self.theorem_symbol != PINNED_DEFECT_THEOREM:
            raise ValueError("theorem_symbol must be requested_cross_defect")
        if self.dependency_symbols != REQUIRED_DEFECT_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned defect theorem chain exactly")

        for name in (
            "canonical_active_family_identity_certified",
            "actual_strip_membership_certified",
            "actual_cross_identity_certified",
            "requested_stress_identity_certified",
            "physical_scale_identity_certified",
            "prepared_N_identity_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class ActualSignedMeanDefectAdmission:
    """Compose canonical signed scope identity with the exact cutoff-defect theorem."""

    canonical_scope: LargeBandCanonicalSignedCommonCurlScopeAdmission
    witness: ActualSignedMeanDefectWitness

    def __post_init__(self) -> None:
        if not isinstance(
            self.canonical_scope, LargeBandCanonicalSignedCommonCurlScopeAdmission
        ):
            raise TypeError(
                "canonical_scope must be a LargeBandCanonicalSignedCommonCurlScopeAdmission"
            )
        if not isinstance(self.witness, ActualSignedMeanDefectWitness):
            raise TypeError("witness must be an ActualSignedMeanDefectWitness")
        if not self.canonical_scope.canonical_to_actual_signed_scope_machine_checked:
            raise ValueError("canonical signed scope is no longer internally coherent")

        export = self.canonical_scope.canonical_export.export
        if (self.witness.B, self.witness.N0) != (export.B, export.N0):
            raise ValueError("mean-defect theorem B/N0 does not match the canonical signed scope")
        if self.witness.prepared_N != export.prepared_N:
            raise ValueError("mean-defect Prepared N does not match the canonical export")
        if self.witness.lean_repository != export.formal_repository:
            raise ValueError("mean-defect repository does not match the canonical export")
        if self.witness.lean_commit != export.formal_commit:
            raise ValueError("mean-defect commit does not match the canonical export")

    @property
    def requested_cross_defect_identity_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def requested_cross_defect_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_mean_cross_values_materialized(self) -> bool:
        return False

    @property
    def missing_weight_values_materialized(self) -> bool:
        return False

    @property
    def finite_head_mean_defect_values_materialized(self) -> bool:
        return False

    @property
    def finite_head_mean_defect_solved(self) -> bool:
        return False

    @property
    def compact_mean_correction_available(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
