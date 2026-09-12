"""Fail-closed admission for the two conclusions of Section 7 ``phase_estimates``.

The large-band path already binds an active signed family index, its carrier/
slot membership, the sign-free base-source revision, and the scalar
``phaseConstant(M) / S_*`` envelope.  What it did not yet type was the claim
that a theorem witness is actually an application of the pinned
``BasePhaseGeometry.FamilyData.phase_estimates`` theorem for that same family.

This module closes only that interface gap.  It admits theorem-provenanced
witness metadata for the theorem's two conclusions on one active carrier/slot:

* the actual phase normal is close to the frozen packed base normal; and
* the actual ``phase.velocity`` is bounded by the same family envelope.

The witness must point to the exact pinned Lean repository/commit/theorem and
to the exact signed label and base-source revision already admitted upstream.
Sampled, fitted, and numeric-scan evidence is rejected.  No vector field is
constructed or sampled here, and Python does not verify the upstream theorem
proof itself.  Passing this gate therefore remains ``formal-structure`` and
cannot promote Eqs. (7.9)--(7.11) or paper-exact velocity truth.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_large_band_family_inputs import LargeBandFamilyInputAdmission
from .phase_large_band_local_base import AsymptoticSlowLabel


PINNED_LEAN_REPOSITORY = "openai/NavierStokesAndEuler"
PINNED_LEAN_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_PHASE_ESTIMATES_THEOREM = "BasePhaseGeometry.FamilyData.phase_estimates"
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
class LargeBandPhaseEstimatesWitness:
    """Theorem-facing identity witness for one ``phase_estimates`` application.

    The four certification flags are intentionally qualitative theorem facts,
    not booleans inferred from numerical samples.  In particular, the normal
    identity and ``phase.velocity`` identity prevent a caller from substituting
    a convenient surrogate vector while keeping the same scalar envelope.
    """

    label: AsymptoticSlowLabel
    source_id: str
    source_revision: str
    lean_repository: str
    lean_commit: str
    theorem_name: str
    evidence_kind: str
    provenance: str
    theorem_application_certified: bool
    actual_normal_identity_certified: bool
    actual_phase_velocity_identity_certified: bool
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
        object.__setattr__(self, "theorem_name", _nonempty_text(self.theorem_name, "theorem_name"))
        object.__setattr__(
            self, "provenance", _nonempty_text(self.provenance, "theorem provenance")
        )
        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; sampled/fitted evidence is rejected"
            )
        if self.lean_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("lean_repository must match the pinned formal source")
        if self.lean_commit != PINNED_LEAN_COMMIT:
            raise ValueError("lean_commit must match the pinned formal source")
        if self.theorem_name != PINNED_PHASE_ESTIMATES_THEOREM:
            raise ValueError("theorem_name must identify the pinned phase_estimates theorem")
        for name in (
            "theorem_application_certified",
            "actual_normal_identity_certified",
            "actual_phase_velocity_identity_certified",
            "uniform_carrier_slot_conclusion_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def theorem_key(self) -> tuple[str, str, str]:
        return self.lean_repository, self.lean_commit, self.theorem_name


@dataclass(frozen=True)
class LargeBandPhaseEstimatesAdmission:
    """Bind a pinned theorem witness to one admitted active family index.

    A successful object exposes the *theorem envelope* for both Lean
    conclusions.  It does not evaluate either vector expression and therefore
    does not claim that the theorem application has been machine-checked in
    this Python repository.
    """

    family: LargeBandFamilyInputAdmission
    witness: LargeBandPhaseEstimatesWitness

    def __post_init__(self) -> None:
        if not isinstance(self.family, LargeBandFamilyInputAdmission):
            raise TypeError("family must be a LargeBandFamilyInputAdmission")
        if not isinstance(self.witness, LargeBandPhaseEstimatesWitness):
            raise TypeError("witness must be a LargeBandPhaseEstimatesWitness")
        failed = [name for name, ok in self.admission_checks().items() if not ok]
        if failed:
            raise ValueError("uncertified large-band phase_estimates admission: " + ", ".join(failed))

    @property
    def label(self) -> AsymptoticSlowLabel:
        return self.family.label

    @property
    def source_key(self) -> tuple[str, str]:
        return self.family.source_key

    @property
    def theorem_envelope(self) -> float:
        """The pinned common bound ``phaseConstant(M) / S_*``.

        The value is inherited from the already-admitted scalar family bridge,
        not measured from any reconstructed field.
        """

        bound = float(self.family.family_phase_bound)
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("family phase theorem envelope must be finite and positive")
        return bound

    @property
    def normal_error_bound(self) -> float:
        return self.theorem_envelope

    @property
    def phase_velocity_bound(self) -> float:
        return self.theorem_envelope

    def admission_checks(self) -> dict[str, bool]:
        w = self.witness
        return {
            "signed_label_identity": w.label == self.family.label,
            "source_revision_identity": w.source_key == self.family.source_key,
            "family_input_admission_intact": all(self.family.admission_checks().values()),
            "pinned_lean_repository": w.lean_repository == PINNED_LEAN_REPOSITORY,
            "pinned_lean_commit": w.lean_commit == PINNED_LEAN_COMMIT,
            "pinned_phase_estimates_theorem": w.theorem_name == PINNED_PHASE_ESTIMATES_THEOREM,
            "theorem_application_certified": w.theorem_application_certified is True,
            "actual_normal_identity_certified": w.actual_normal_identity_certified is True,
            "actual_phase_velocity_identity_certified": w.actual_phase_velocity_identity_certified is True,
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
    def actual_phase_vectors_materialized(self) -> bool:
        return False

    @property
    def family_phase_estimates_machine_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
