"""Fail-closed Section 7 bridge for pinned ``base_estimates`` -> ``coordinate_errors``.

The large-band path already admits one active signed family, binds it to a
stable sign-free base-source revision, and types an application of the pinned
``BasePhaseGeometry.FamilyData.phase_estimates`` theorem.  The next Lean step
uses ``FamilyData.base_estimates`` together with those phase estimates to prove
``FamilyData.coordinate_errors`` for the moving-frame coefficients.

This module records only that theorem/application identity.  A witness must
name the exact pinned Lean repository, commit, and both theorem names; reuse the
exact signed label and base-source revision already admitted upstream; and
certify the actual base/shear/frame expressions rather than substitute sampled
or fitted surrogates.

No base vector, shear vector, moving-frame coefficient, partition, or wave is
materialized here.  The scalar bounds exposed below are theorem envelopes
inherited from the already-landed large-band geometry.  Passing this gate is
therefore only ``formal-structure`` and does not establish Eqs. (7.9)-(7.11)
for the manuscript construction.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_phase_estimates import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
    LargeBandPhaseEstimatesAdmission,
)


PINNED_BASE_ESTIMATES_THEOREM = "BasePhaseGeometry.FamilyData.base_estimates"
PINNED_COORDINATE_ERRORS_THEOREM = "BasePhaseGeometry.FamilyData.coordinate_errors"
_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class LargeBandCoordinateErrorsWitness:
    """Theorem-facing identity witness for ``base_estimates``/``coordinate_errors``.

    The certification flags are qualitative theorem facts, not observations.
    They force the downstream bridge to refer to the actual base value, actual
    shear-vector expression, and actual moving-frame ``errorA/B/C`` expressions
    appearing in the pinned formal source.
    """

    label: AsymptoticSlowLabel
    source_id: str
    source_revision: str
    lean_repository: str
    lean_commit: str
    base_estimates_theorem_name: str
    coordinate_errors_theorem_name: str
    evidence_kind: str
    provenance: str
    base_estimates_application_certified: bool
    actual_base_value_identity_certified: bool
    actual_shear_vector_identity_certified: bool
    uniform_carrier_base_conclusion_certified: bool
    coordinate_errors_application_certified: bool
    actual_frame_errors_identity_certified: bool
    phase_estimates_dependency_certified: bool
    base_estimates_dependency_certified: bool
    uniform_carrier_slot_conclusion_certified: bool

    def __post_init__(self) -> None:
        if not isinstance(self.label, AsymptoticSlowLabel):
            raise TypeError("label must be an AsymptoticSlowLabel")
        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self, "source_revision", _nonempty_text(self.source_revision, "source_revision")
        )
        object.__setattr__(
            self, "lean_repository", _nonempty_text(self.lean_repository, "lean_repository")
        )
        object.__setattr__(self, "lean_commit", _nonempty_text(self.lean_commit, "lean_commit"))
        object.__setattr__(
            self,
            "base_estimates_theorem_name",
            _nonempty_text(self.base_estimates_theorem_name, "base_estimates_theorem_name"),
        )
        object.__setattr__(
            self,
            "coordinate_errors_theorem_name",
            _nonempty_text(self.coordinate_errors_theorem_name, "coordinate_errors_theorem_name"),
        )
        object.__setattr__(
            self, "provenance", _nonempty_text(self.provenance, "theorem provenance")
        )

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted evidence is rejected"
            )
        if self.lean_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("lean_repository must match the pinned formal source")
        if self.lean_commit != PINNED_LEAN_COMMIT:
            raise ValueError("lean_commit must match the pinned formal source")
        if self.base_estimates_theorem_name != PINNED_BASE_ESTIMATES_THEOREM:
            raise ValueError(
                "base_estimates_theorem_name must identify the pinned base_estimates theorem"
            )
        if self.coordinate_errors_theorem_name != PINNED_COORDINATE_ERRORS_THEOREM:
            raise ValueError(
                "coordinate_errors_theorem_name must identify the pinned coordinate_errors theorem"
            )

        for name in (
            "base_estimates_application_certified",
            "actual_base_value_identity_certified",
            "actual_shear_vector_identity_certified",
            "uniform_carrier_base_conclusion_certified",
            "coordinate_errors_application_certified",
            "actual_frame_errors_identity_certified",
            "phase_estimates_dependency_certified",
            "base_estimates_dependency_certified",
            "uniform_carrier_slot_conclusion_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def theorem_key(self) -> tuple[str, str, str, str]:
        return (
            self.lean_repository,
            self.lean_commit,
            self.base_estimates_theorem_name,
            self.coordinate_errors_theorem_name,
        )


@dataclass(frozen=True)
class LargeBandCoordinateErrorsAdmission:
    """Bind the two pinned downstream theorem applications to one phase admission."""

    phase: LargeBandPhaseEstimatesAdmission
    witness: LargeBandCoordinateErrorsWitness

    def __post_init__(self) -> None:
        if not isinstance(self.phase, LargeBandPhaseEstimatesAdmission):
            raise TypeError("phase must be a LargeBandPhaseEstimatesAdmission")
        if not isinstance(self.witness, LargeBandCoordinateErrorsWitness):
            raise TypeError("witness must be a LargeBandCoordinateErrorsWitness")
        failed = [name for name, ok in self.admission_checks().items() if not ok]
        if failed:
            raise ValueError(
                "uncertified large-band base_estimates/coordinate_errors admission: "
                + ", ".join(failed)
            )

    @property
    def label(self) -> AsymptoticSlowLabel:
        return self.phase.label

    @property
    def source_key(self) -> tuple[str, str]:
        return self.phase.source_key

    @property
    def base_theorem_envelope(self) -> float:
        """Pinned ``16*M^2 / D.scale(i)`` envelope with ``D.scale(i)=S_*`` here."""

        M = float(self.phase.family.witness.M)
        S_star = float(self.phase.family.frame.S_star)
        bound = 16.0 * M * M / S_star
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("base_estimates theorem envelope must be finite and positive")
        return bound

    @property
    def base_value_error_bound(self) -> float:
        return self.base_theorem_envelope

    @property
    def base_shear_vector_error_bound(self) -> float:
        return self.base_theorem_envelope

    @property
    def coordinate_theorem_envelope(self) -> float:
        """Pinned ``coordinateConstant(M,u) / D.scale(i)`` scalar envelope."""

        bound = float(self.phase.family.frame.frame_error_bound)
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("coordinate_errors theorem envelope must be finite and positive")
        return bound

    @property
    def coordinate_error_A_bound(self) -> float:
        return self.coordinate_theorem_envelope

    @property
    def coordinate_error_B_bound(self) -> float:
        return self.coordinate_theorem_envelope

    @property
    def coordinate_error_C_bound(self) -> float:
        return self.coordinate_theorem_envelope

    def admission_checks(self) -> dict[str, bool]:
        w = self.witness
        return {
            "signed_label_identity": w.label == self.phase.label,
            "source_revision_identity": w.source_key == self.phase.source_key,
            "phase_estimates_admission_intact": all(self.phase.admission_checks().values()),
            "pinned_lean_repository": w.lean_repository == PINNED_LEAN_REPOSITORY,
            "pinned_lean_commit": w.lean_commit == PINNED_LEAN_COMMIT,
            "pinned_base_estimates_theorem": (
                w.base_estimates_theorem_name == PINNED_BASE_ESTIMATES_THEOREM
            ),
            "pinned_coordinate_errors_theorem": (
                w.coordinate_errors_theorem_name == PINNED_COORDINATE_ERRORS_THEOREM
            ),
            "base_estimates_application_certified": (
                w.base_estimates_application_certified is True
            ),
            "actual_base_value_identity_certified": (
                w.actual_base_value_identity_certified is True
            ),
            "actual_shear_vector_identity_certified": (
                w.actual_shear_vector_identity_certified is True
            ),
            "uniform_carrier_base_conclusion_certified": (
                w.uniform_carrier_base_conclusion_certified is True
            ),
            "coordinate_errors_application_certified": (
                w.coordinate_errors_application_certified is True
            ),
            "actual_frame_errors_identity_certified": (
                w.actual_frame_errors_identity_certified is True
            ),
            "phase_estimates_dependency_certified": (
                w.phase_estimates_dependency_certified is True
            ),
            "base_estimates_dependency_certified": (
                w.base_estimates_dependency_certified is True
            ),
            "uniform_carrier_slot_conclusion_certified": (
                w.uniform_carrier_slot_conclusion_certified is True
            ),
            "theorem_evidence_not_sampled": w.evidence_kind in _ACCEPTED_EVIDENCE,
            "provenance_present": bool(w.provenance),
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_base_values_materialized(self) -> bool:
        return False

    @property
    def actual_shear_vectors_materialized(self) -> bool:
        return False

    @property
    def actual_coordinate_errors_materialized(self) -> bool:
        return False

    @property
    def base_estimates_machine_verified(self) -> bool:
        return False

    @property
    def coordinate_errors_machine_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
