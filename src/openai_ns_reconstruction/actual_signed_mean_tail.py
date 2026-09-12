"""Bind the canonical signed scope to the pinned exact mean-cross tail theorem.

Section 8 mean repair must eventually consume the *actual* signed correction,
not a fitted covariance surrogate.  The pinned formal source already proves a
strong prerequisite: ``ActualSignedMeanBinding.literal_requested_cross_tail``
identifies the mean cross of the literal correction-cycle signed field with the
requested stress on every tail band satisfying ``Prepared.N + 1 <= n``.

This module adds only a fail-closed admission boundary for that theorem
application.  It composes the already-landed canonical signed-scope identity
with one externally certified theorem application, requiring the same ``B/N0``
and selected canonical ``Prepared.N``.  Lean cycle states and points remain
opaque strings; Python does not encode or reconstruct them.

Successful admission is still ``formal-structure``.  No Lean proof term is
replayed here, no mean field is materialized, no finite-head defect is solved,
and no compact mean correction or paper-exact velocity is promoted.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .actual_signed_canonical_scope import (
    LargeBandCanonicalSignedCommonCurlScopeAdmission,
)
from .actual_signed_common_curl import PINNED_LEAN_COMMIT, PINNED_LEAN_REPOSITORY


PINNED_MEAN_FILE = "NavierStokes/ActualSignedMeanBinding.lean"
PINNED_LITERAL_TAIL_THEOREM = (
    "NavierStokes.ActualSignedMeanBinding.literal_requested_cross_tail"
)
PINNED_REQUESTED_TAIL_THEOREM = (
    "NavierStokes.ActualSignedMeanBinding.requested_cross_tail"
)
PINNED_CYCLE_CROSS_THEOREM = "NavierStokes.ActualSignedMeanBinding.cycle_cross_eq"
PINNED_PARAMETERS_EQ_CYCLE = (
    "NavierStokes.ActualSignedMeanBinding.parameters_eq_cycle"
)
REQUIRED_MEAN_SYMBOLS = (
    PINNED_LITERAL_TAIL_THEOREM,
    PINNED_REQUESTED_TAIL_THEOREM,
    PINNED_CYCLE_CROSS_THEOREM,
    PINNED_PARAMETERS_EQ_CYCLE,
)


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
class ActualSignedMeanTailWitness:
    """Opaque metadata for one literal ``requested_cross_tail`` application."""

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
    tail_band_hypothesis_certified: bool
    actual_strip_membership_certified: bool
    literal_cycle_parameters_identity_certified: bool
    theorem_application_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_MEAN_FILE
    theorem_symbol: str = PINNED_LITERAL_TAIL_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_MEAN_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "prepared_N", _natural(self.prepared_N, "prepared_N"))
        object.__setattr__(self, "band", _natural(self.band, "band"))
        component = _natural(self.component, "component")
        if component not in (0, 1):
            raise ValueError("component must be 0 (theta) or 1 (axial)")
        object.__setattr__(self, "component", component)
        object.__setattr__(self, "cycle_state_repr", _nonempty_text(self.cycle_state_repr, "cycle_state_repr"))
        object.__setattr__(self, "point_repr", _nonempty_text(self.point_repr, "point_repr"))
        object.__setattr__(self, "application_id", _nonempty_text(self.application_id, "application_id"))
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
        if self.theorem_symbol != PINNED_LITERAL_TAIL_THEOREM:
            raise ValueError("theorem_symbol must be literal_requested_cross_tail")
        if self.dependency_symbols != REQUIRED_MEAN_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned theorem chain exactly")

        if self.band < self.prepared_N + 1:
            raise ValueError("band must satisfy the pinned tail hypothesis Prepared.N + 1 <= n")
        for name in (
            "active_labels_identity_certified",
            "tail_band_hypothesis_certified",
            "actual_strip_membership_certified",
            "literal_cycle_parameters_identity_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class ActualSignedMeanTailAdmission:
    """Compose canonical signed scope identity with exact tail mean cancellation."""

    canonical_scope: LargeBandCanonicalSignedCommonCurlScopeAdmission
    witness: ActualSignedMeanTailWitness

    def __post_init__(self) -> None:
        if not isinstance(
            self.canonical_scope, LargeBandCanonicalSignedCommonCurlScopeAdmission
        ):
            raise TypeError(
                "canonical_scope must be a LargeBandCanonicalSignedCommonCurlScopeAdmission"
            )
        if not isinstance(self.witness, ActualSignedMeanTailWitness):
            raise TypeError("witness must be an ActualSignedMeanTailWitness")
        if not self.canonical_scope.canonical_to_actual_signed_scope_machine_checked:
            raise ValueError("canonical signed scope is no longer internally coherent")

        export = self.canonical_scope.canonical_export.export
        if (self.witness.B, self.witness.N0) != (export.B, export.N0):
            raise ValueError("mean-tail theorem B/N0 does not match the canonical signed scope")
        if self.witness.prepared_N != export.prepared_N:
            raise ValueError("mean-tail Prepared N does not match the canonical export")
        if self.witness.lean_repository != export.formal_repository:
            raise ValueError("mean-tail repository does not match the canonical export")
        if self.witness.lean_commit != export.formal_commit:
            raise ValueError("mean-tail commit does not match the canonical export")

    @property
    def requested_cross_tail_identity_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def literal_requested_cross_tail_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_mean_cross_values_materialized(self) -> bool:
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
