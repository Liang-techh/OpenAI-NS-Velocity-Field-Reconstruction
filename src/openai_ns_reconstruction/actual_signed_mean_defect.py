"""Admit the pinned finite-head signed mean-cross defect identity.

Section 8 cannot replace the low-band mean deficit by a synthetic covariance
surrogate.  The pinned formal source proves the exact all-band identity
``ActualSignedMeanBinding.requested_cross_defect``:

``meanBar(actualCross) - requestedStress
  = -missingWeight Prepared.N (physicalScale n x) * requestedStress``.

This module deliberately narrows that identity to the finite head
``n <= Prepared.N``.  The complementary tail was admitted separately by
:mod:`actual_signed_mean_tail`, where exact cancellation is certified for
``Prepared.N + 1 <= n``.

Only theorem-grade metadata is accepted.  The requested stress, missing-weight
factor, cycle state, and point remain opaque; Python neither samples nor
materializes the debt and therefore cannot yet feed it into the compact mean
repair as paper-exact data.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .actual_signed_canonical_scope import (
    LargeBandCanonicalSignedCommonCurlScopeAdmission,
)
from .actual_signed_common_curl import PINNED_LEAN_COMMIT, PINNED_LEAN_REPOSITORY
from .actual_signed_mean_tail import PINNED_MEAN_FILE


PINNED_COVARIANCE_FILE = "NavierStokes/ActualPrimaryCovariance.lean"
PINNED_REQUESTED_CROSS_DEFECT_THEOREM = (
    "NavierStokes.ActualSignedMeanBinding.requested_cross_defect"
)
PINNED_REQUESTED_CROSS_FACTOR_THEOREM = (
    "NavierStokes.ActualSignedMeanBinding.requested_cross_factor"
)
PINNED_PARTITION_DEFECT_THEOREM = (
    "NavierStokes.ActualPrimaryCovariance.partitionFactor_eq_one_sub_missing"
)
REQUIRED_DEFECT_SYMBOLS = (
    PINNED_REQUESTED_CROSS_DEFECT_THEOREM,
    PINNED_REQUESTED_CROSS_FACTOR_THEOREM,
    PINNED_PARTITION_DEFECT_THEOREM,
)
REQUIRED_DEFECT_FILES = (PINNED_MEAN_FILE, PINNED_COVARIANCE_FILE)


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
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class ActualSignedMeanDefectWitness:
    """Opaque metadata for one finite-head ``requested_cross_defect`` application."""

    B: int
    N0: int
    prepared_N: int
    band: int
    component: int
    cycle_state_repr: str
    point_repr: str
    application_id: str
    evidence_kind: str
    provenance: str
    active_labels_identity_certified: bool
    finite_head_band_certified: bool
    actual_strip_membership_certified: bool
    requested_stress_identity_certified: bool
    missing_weight_factor_identity_certified: bool
    theorem_application_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_MEAN_FILE
    theorem_symbol: str = PINNED_REQUESTED_CROSS_DEFECT_THEOREM
    dependency_files: tuple[str, ...] = REQUIRED_DEFECT_FILES
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
        object.__setattr__(
            self, "cycle_state_repr", _nonempty_text(self.cycle_state_repr, "cycle_state_repr")
        )
        object.__setattr__(self, "point_repr", _nonempty_text(self.point_repr, "point_repr"))
        object.__setattr__(
            self, "application_id", _nonempty_text(self.application_id, "application_id")
        )
        object.__setattr__(self, "provenance", _nonempty_text(self.provenance, "provenance"))

        if self.evidence_kind != "formal-theorem":
            raise ValueError(
                "evidence_kind must be formal-theorem; sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.lean_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("lean_repository must match the pinned formal source exactly")
        if self.lean_commit != PINNED_LEAN_COMMIT:
            raise ValueError("lean_commit must match the pinned formal source exactly")
        if self.lean_file != PINNED_MEAN_FILE:
            raise ValueError("lean_file must match the pinned actual signed mean source exactly")
        if self.theorem_symbol != PINNED_REQUESTED_CROSS_DEFECT_THEOREM:
            raise ValueError("theorem_symbol must be requested_cross_defect")
        if self.dependency_files != REQUIRED_DEFECT_FILES:
            raise ValueError("dependency_files must match the pinned defect source chain exactly")
        if self.dependency_symbols != REQUIRED_DEFECT_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned defect theorem chain exactly")

        # requested_cross_defect itself is all-band.  We intentionally admit it
        # here only on the finite complement of the previously admitted tail.
        if self.band > self.prepared_N:
            raise ValueError("finite-head admission requires n <= Prepared.N")
        for name in (
            "active_labels_identity_certified",
            "finite_head_band_certified",
            "actual_strip_membership_certified",
            "requested_stress_identity_certified",
            "missing_weight_factor_identity_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class ActualSignedMeanDefectAdmission:
    """Bind the exact finite-head defect identity to the canonical signed scope."""

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
    def finite_head_band_admitted(self) -> bool:
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
    def finite_head_mean_debt_materialized(self) -> bool:
        return False

    @property
    def compact_mean_correction_available(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
