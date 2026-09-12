"""Fail-closed bridge from large-band base sources to the formal Prepared constructor.

The large-band Section 6/7 adapters already bind each sign-free slow box to a
stable ``(source_id, source_revision)``.  That identity alone does not say that
the source is the one constructed by the pinned formal primary-geometry
assembly.

The pinned ``PrimaryGeometryAssembly.lean`` source supplies a stronger chain:
``Prepared`` contains the LocalBase family and polynomial jets;
``exists_prepared`` constructs such data from the actual slow base and true
cone; ``prepared`` selects that witness; ``family`` reuses ``a.base`` and only
then inserts ``phaseSign c``; and ``phases`` builds both ``Fin 2`` branches from
the same selected ``prepared`` object.

This module records a theorem-provenanced identity bridge to that chain.  It is
intentionally metadata-only: Python does not evaluate Lean's noncomputable
objects, discharge ``exists_prepared`` hypotheses, or materialize any base or
phase field.  Sampled, fitted, and numeric-scan evidence is rejected.  Until an
actual exported theorem witness is connected, the status remains strictly
``formal-structure`` and ``paper_exact_velocity_available`` is false.
"""
from __future__ import annotations

from dataclasses import dataclass

from .phase_large_band_base_source import LargeBandBaseSourceBinding


PINNED_LEAN_REPOSITORY = "openai/NavierStokesAndEuler"
PINNED_LEAN_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_LEAN_FILE = "NavierStokes/PrimaryGeometryAssembly.lean"
PINNED_PREPARED_DECL = "PrimaryGeometryAssembly.Prepared"
PINNED_EXISTS_PREPARED = "PrimaryGeometryAssembly.exists_prepared"
PINNED_PREPARED_DEF = "PrimaryGeometryAssembly.prepared"
PINNED_FAMILY_DEF = "PrimaryGeometryAssembly.family"
PINNED_CONSTRUCTION_DEF = "PrimaryGeometryAssembly.construction"
PINNED_PHASES_DEF = "PrimaryGeometryAssembly.phases"

_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _pinned(value: object, expected: str, name: str) -> str:
    text = _nonempty_text(value, name)
    if text != expected:
        raise ValueError(f"{name} must match the pinned formal source exactly")
    return text


@dataclass(frozen=True)
class PrimaryGeometryPreparedSourceWitness:
    """External theorem identity for one selected formal ``Prepared`` source.

    ``prepared_instance_id`` is an opaque stable identifier supplied by the
    theorem/export layer.  It distinguishes concrete applications of the same
    pinned constructor without pretending that Python can evaluate the
    noncomputable Lean object.

    The five certification flags are deliberately separate.  In particular,
    existence of a ``Prepared`` object is not enough: downstream use must also
    certify that the selected ``prepared`` is the same source exposed through
    ``family ... .base`` and reused by the two ``phases`` branches.
    """

    source_id: str
    source_revision: str
    prepared_instance_id: str
    evidence_kind: str
    provenance: str
    exists_prepared_application_certified: bool
    prepared_choice_identity_certified: bool
    family_base_source_identity_certified: bool
    phases_reuse_prepared_certified: bool
    both_signs_share_prepared_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_LEAN_FILE
    prepared_decl: str = PINNED_PREPARED_DECL
    exists_prepared_theorem: str = PINNED_EXISTS_PREPARED
    prepared_def: str = PINNED_PREPARED_DEF
    family_def: str = PINNED_FAMILY_DEF
    construction_def: str = PINNED_CONSTRUCTION_DEF
    phases_def: str = PINNED_PHASES_DEF

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self, "source_revision", _nonempty_text(self.source_revision, "source_revision")
        )
        object.__setattr__(
            self,
            "prepared_instance_id",
            _nonempty_text(self.prepared_instance_id, "prepared_instance_id"),
        )
        object.__setattr__(
            self, "provenance", _nonempty_text(self.provenance, "theorem provenance")
        )
        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted evidence is rejected"
            )
        for name in (
            "exists_prepared_application_certified",
            "prepared_choice_identity_certified",
            "family_base_source_identity_certified",
            "phases_reuse_prepared_certified",
            "both_signs_share_prepared_certified",
        ):
            _strict_true(getattr(self, name), name)

        pins = (
            ("lean_repository", PINNED_LEAN_REPOSITORY),
            ("lean_commit", PINNED_LEAN_COMMIT),
            ("lean_file", PINNED_LEAN_FILE),
            ("prepared_decl", PINNED_PREPARED_DECL),
            ("exists_prepared_theorem", PINNED_EXISTS_PREPARED),
            ("prepared_def", PINNED_PREPARED_DEF),
            ("family_def", PINNED_FAMILY_DEF),
            ("construction_def", PINNED_CONSTRUCTION_DEF),
            ("phases_def", PINNED_PHASES_DEF),
        )
        for name, expected in pins:
            object.__setattr__(self, name, _pinned(getattr(self, name), expected, name))

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.source_id, self.source_revision, self.prepared_instance_id


@dataclass(frozen=True)
class LargeBandPreparedSourceBinding:
    """Bind an admitted sign-free box source to the pinned ``Prepared`` chain."""

    base_source: LargeBandBaseSourceBinding
    prepared: PrimaryGeometryPreparedSourceWitness

    def __post_init__(self) -> None:
        if not isinstance(self.base_source, LargeBandBaseSourceBinding):
            raise TypeError("base_source must be a LargeBandBaseSourceBinding")
        if not isinstance(self.prepared, PrimaryGeometryPreparedSourceWitness):
            raise TypeError("prepared must be a PrimaryGeometryPreparedSourceWitness")
        if self.prepared.source_key != self.base_source.source_key:
            raise ValueError(
                "Prepared source identity must exactly match the admitted base-source revision"
            )
        failed = [name for name, ok in self.binding_checks().items() if not ok]
        if failed:
            raise ValueError(
                "uncertified PrimaryGeometry Prepared source binding: " + ", ".join(failed)
            )

    @property
    def source_key(self) -> tuple[str, str]:
        return self.base_source.source_key

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.prepared.prepared_key

    @property
    def sign_pair(self):
        """Return the exact two signed labels already tied to the common source."""

        return self.base_source.sign_pair

    def binding_checks(self) -> dict[str, bool]:
        return {
            "base_source_binding_intact": all(self.base_source.binding_checks().values()),
            "source_revision_identity": self.prepared.source_key == self.base_source.source_key,
            "exists_prepared_application_certified": (
                self.prepared.exists_prepared_application_certified is True
            ),
            "prepared_choice_identity_certified": (
                self.prepared.prepared_choice_identity_certified is True
            ),
            "family_base_source_identity_certified": (
                self.prepared.family_base_source_identity_certified is True
            ),
            "phases_reuse_prepared_certified": self.prepared.phases_reuse_prepared_certified is True,
            "both_signs_share_prepared_certified": (
                self.prepared.both_signs_share_prepared_certified is True
            ),
            "theorem_evidence_not_sampled": self.prepared.evidence_kind in _ACCEPTED_EVIDENCE,
            "prepared_instance_identity_present": bool(self.prepared.prepared_instance_id),
            "provenance_present": bool(self.prepared.provenance),
            "formal_source_pin_intact": (
                self.prepared.lean_repository == PINNED_LEAN_REPOSITORY
                and self.prepared.lean_commit == PINNED_LEAN_COMMIT
                and self.prepared.lean_file == PINNED_LEAN_FILE
                and self.prepared.prepared_decl == PINNED_PREPARED_DECL
                and self.prepared.exists_prepared_theorem == PINNED_EXISTS_PREPARED
                and self.prepared.prepared_def == PINNED_PREPARED_DEF
                and self.prepared.family_def == PINNED_FAMILY_DEF
                and self.prepared.construction_def == PINNED_CONSTRUCTION_DEF
                and self.prepared.phases_def == PINNED_PHASES_DEF
            ),
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_prepared_application_machine_verified(self) -> bool:
        return False

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
