"""Fail-closed Section 7 admission for the pinned ``FamilyData.damping_error`` theorem.

The large-band path already admits one active signed family index, binds it to
one stable base-source revision, and checks the scalar frame/damping
specialization without forming the underflow-prone chart scale ``Q=2^-ell``.
The pinned Lean theorem

``BasePhaseGeometry.FamilyData.damping_error``

then identifies the *actual* viscosity-damping discrepancy on an active carrier
with the uniform envelope ``dampingConstant(M) / D.scale(i)``.

This module records only that theorem/application identity.  It does not
materialize the damping field or infer the theorem from sampled values.  A
witness must name the exact pinned Lean source, reuse the already-admitted
signed family/source identity, and theorem-certify the actual and reference
damping expressions, the ``D.scale(i)=S_*`` identity, and the reference
viscosity/radius identities.

Passing this gate is therefore only ``formal-structure``.  In particular it
does not machine-prove the manuscript family, Eqs. (7.9)-(7.11), or a
paper-exact velocity field.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_frame_bounds import damping_constant
from .phase_large_band_family_inputs import LargeBandFamilyInputAdmission
from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_phase_estimates import PINNED_LEAN_COMMIT, PINNED_LEAN_REPOSITORY


PINNED_DAMPING_ERROR_THEOREM = "BasePhaseGeometry.FamilyData.damping_error"
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
class LargeBandDampingErrorWitness:
    """Theorem-facing identity witness for ``FamilyData.damping_error``.

    Every boolean below is a theorem/provenance fact, never a sampled
    observation.  The witness deliberately carries no damping samples or
    user-adjustable error tolerance.
    """

    label: AsymptoticSlowLabel
    source_id: str
    source_revision: str
    lean_repository: str
    lean_commit: str
    theorem_name: str
    evidence_kind: str
    provenance: str
    damping_error_application_certified: bool
    carrier_membership_certified: bool
    actual_damping_term_identity_certified: bool
    reference_damping_term_identity_certified: bool
    viscosity_damping_family_identity_certified: bool
    reference_viscosity_identity_certified: bool
    reference_radius_center_identity_certified: bool
    scale_identity_certified: bool
    family_input_dependency_certified: bool

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
        object.__setattr__(self, "theorem_name", _nonempty_text(self.theorem_name, "theorem_name"))
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
        if self.theorem_name != PINNED_DAMPING_ERROR_THEOREM:
            raise ValueError("theorem_name must identify the pinned damping_error theorem")

        for name in (
            "damping_error_application_certified",
            "carrier_membership_certified",
            "actual_damping_term_identity_certified",
            "reference_damping_term_identity_certified",
            "viscosity_damping_family_identity_certified",
            "reference_viscosity_identity_certified",
            "reference_radius_center_identity_certified",
            "scale_identity_certified",
            "family_input_dependency_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def theorem_key(self) -> tuple[str, str, str]:
        return self.lean_repository, self.lean_commit, self.theorem_name


@dataclass(frozen=True)
class LargeBandDampingErrorAdmission:
    """Bind the pinned damping theorem application to one admitted family."""

    family: LargeBandFamilyInputAdmission
    witness: LargeBandDampingErrorWitness

    def __post_init__(self) -> None:
        if not isinstance(self.family, LargeBandFamilyInputAdmission):
            raise TypeError("family must be a LargeBandFamilyInputAdmission")
        if not isinstance(self.witness, LargeBandDampingErrorWitness):
            raise TypeError("witness must be a LargeBandDampingErrorWitness")
        failed = [name for name, ok in self.admission_checks().items() if not ok]
        if failed:
            raise ValueError(
                "uncertified large-band damping_error admission: " + ", ".join(failed)
            )

    @property
    def label(self) -> AsymptoticSlowLabel:
        return self.family.label

    @property
    def source_key(self) -> tuple[str, str]:
        return self.family.source_key

    @property
    def damping_theorem_envelope(self) -> float:
        """Pinned ``dampingConstant(M) / D.scale(i)`` with ``D.scale(i)=S_*``."""

        bound = float(self.family.frame.damping_error_bound)
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("damping_error theorem envelope must be finite and positive")
        return bound

    @property
    def constant_over_scale_envelope(self) -> float:
        """Independent named-constant evaluation of the pinned theorem envelope."""

        M = float(self.family.witness.M)
        S_star = float(self.family.frame.S_star)
        bound = damping_constant(M) / S_star
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("dampingConstant(M)/S_* must be finite and positive")
        return bound

    def admission_checks(self) -> dict[str, bool]:
        w = self.witness
        named_constant_matches = math.isclose(
            self.family.frame.damping_error_bound,
            self.constant_over_scale_envelope,
            rel_tol=1e-12,
            abs_tol=0.0,
        )
        return {
            "signed_label_identity": w.label == self.family.label,
            "source_revision_identity": w.source_key == self.family.source_key,
            "family_input_admission_intact": all(self.family.admission_checks().values()),
            "pinned_lean_repository": w.lean_repository == PINNED_LEAN_REPOSITORY,
            "pinned_lean_commit": w.lean_commit == PINNED_LEAN_COMMIT,
            "pinned_damping_error_theorem": w.theorem_name == PINNED_DAMPING_ERROR_THEOREM,
            "damping_error_application_certified": w.damping_error_application_certified is True,
            "carrier_membership_certified": w.carrier_membership_certified is True,
            "actual_damping_term_identity_certified": (
                w.actual_damping_term_identity_certified is True
            ),
            "reference_damping_term_identity_certified": (
                w.reference_damping_term_identity_certified is True
            ),
            "viscosity_damping_family_identity_certified": (
                w.viscosity_damping_family_identity_certified is True
            ),
            "reference_viscosity_identity_certified": (
                w.reference_viscosity_identity_certified is True
            ),
            "reference_radius_center_identity_certified": (
                w.reference_radius_center_identity_certified is True
            ),
            "scale_identity_certified": w.scale_identity_certified is True,
            "family_input_dependency_certified": w.family_input_dependency_certified is True,
            "named_damping_constant_over_scale_identity": named_constant_matches,
            "theorem_evidence_not_sampled": w.evidence_kind in _ACCEPTED_EVIDENCE,
            "provenance_present": bool(w.provenance),
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_damping_field_materialized(self) -> bool:
        return False

    @property
    def damping_error_machine_verified(self) -> bool:
        return False

    @property
    def all_active_family_damping_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
