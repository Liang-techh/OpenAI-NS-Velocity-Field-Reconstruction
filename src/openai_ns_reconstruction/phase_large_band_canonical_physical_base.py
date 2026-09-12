"""Bind a Section 7 physical base to the pinned canonical Prepared construction.

The generic large-band ``Prepared`` bridge intentionally accepts any theorem-
certified application of ``PrimaryGeometryAssembly.prepared``.  The pinned
formal source also exposes the construction actually intended for the primary
geometry:

``canonicalPrepared`` fixes ``upper = 2 * NominalConeAssembly.activeRight W``
and ``canonicalPhases`` builds both ``Fin 2`` branches from that same selected
object.  Since ``frequency`` and ``axial`` use ``FinalSlowBase.scales H v upper
B``, this canonical choice also removes an otherwise unconstrained ``upper``
from the physical-base identity chain.

This module composes that canonical identity with an already admitted
``Prepared -> FinalSlowBase.velocity`` binding.  It is deliberately fail-closed
metadata: Python does not evaluate the noncomputable Lean object, construct the
``FullTrueCone`` witness, or replay the cited definitions.  Sampled, fitted and
numeric-scan evidence is rejected, and successful admission remains
``formal-structure`` until a machine-linked formal export supplies the actual
application.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral

from .phase_large_band_physical_base import LargeBandPhysicalBaseBinding
from .phase_large_band_prepared_source import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_FILE,
    PINNED_LEAN_REPOSITORY,
)


PINNED_CANONICAL_PREPARED_DEF = "PrimaryGeometryAssembly.canonicalPrepared"
PINNED_CANONICAL_PHASES_DEF = "PrimaryGeometryAssembly.canonicalPhases"
PINNED_ACTIVE_RIGHT_DEF = "NominalConeAssembly.activeRight"
PINNED_FREQUENCY_DEF = "PrimaryGeometryAssembly.frequency"
PINNED_AXIAL_DEF = "PrimaryGeometryAssembly.axial"
PINNED_FINAL_SLOW_BASE_SCALES_DEF = "FinalSlowBase.scales"
PINNED_FULL_TRUE_CONE = "LeadingStressWeights.FullTrueCone"

_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _positive_finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite and positive")
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return out


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
class PrimaryGeometryCanonicalPreparedWitness:
    """Theorem-facing identity for one canonical primary-geometry application.

    ``prepared_instance_id`` must identify the same selected object already used
    by the generic Prepared/physical-base chain.  ``B``, ``r0`` and ``N0`` are
    recorded because they are explicit arguments of ``canonicalPrepared`` /
    ``canonicalPhases``.  They are application metadata, not values inferred by
    this module.
    """

    source_id: str
    source_revision: str
    prepared_instance_id: str
    B: int
    r0: float
    N0: int
    evidence_kind: str
    provenance: str
    full_true_cone_witness_certified: bool
    canonical_prepared_application_certified: bool
    canonical_phases_application_certified: bool
    canonical_upper_twice_active_right_certified: bool
    frequency_uses_canonical_scales_certified: bool
    axial_uses_canonical_scales_certified: bool
    both_signs_share_canonical_phases_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_LEAN_FILE
    canonical_prepared_def: str = PINNED_CANONICAL_PREPARED_DEF
    canonical_phases_def: str = PINNED_CANONICAL_PHASES_DEF
    active_right_def: str = PINNED_ACTIVE_RIGHT_DEF
    frequency_def: str = PINNED_FREQUENCY_DEF
    axial_def: str = PINNED_AXIAL_DEF
    final_slow_base_scales_def: str = PINNED_FINAL_SLOW_BASE_SCALES_DEF
    full_true_cone_decl: str = PINNED_FULL_TRUE_CONE

    def __post_init__(self) -> None:
        for name in ("source_id", "source_revision", "prepared_instance_id", "provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "r0", _positive_finite(self.r0, "r0"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        for name in (
            "full_true_cone_witness_certified",
            "canonical_prepared_application_certified",
            "canonical_phases_application_certified",
            "canonical_upper_twice_active_right_certified",
            "frequency_uses_canonical_scales_certified",
            "axial_uses_canonical_scales_certified",
            "both_signs_share_canonical_phases_certified",
        ):
            _strict_true(getattr(self, name), name)

        pins = (
            ("lean_repository", PINNED_LEAN_REPOSITORY),
            ("lean_commit", PINNED_LEAN_COMMIT),
            ("lean_file", PINNED_LEAN_FILE),
            ("canonical_prepared_def", PINNED_CANONICAL_PREPARED_DEF),
            ("canonical_phases_def", PINNED_CANONICAL_PHASES_DEF),
            ("active_right_def", PINNED_ACTIVE_RIGHT_DEF),
            ("frequency_def", PINNED_FREQUENCY_DEF),
            ("axial_def", PINNED_AXIAL_DEF),
            ("final_slow_base_scales_def", PINNED_FINAL_SLOW_BASE_SCALES_DEF),
            ("full_true_cone_decl", PINNED_FULL_TRUE_CONE),
        )
        for name, expected in pins:
            object.__setattr__(self, name, _pinned(getattr(self, name), expected, name))

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.source_id, self.source_revision, self.prepared_instance_id

    @property
    def canonical_parameter_key(self) -> tuple[int, float, int]:
        return self.B, self.r0, self.N0


@dataclass(frozen=True)
class LargeBandCanonicalPhysicalBaseBinding:
    """Require one physical-base binding to use the pinned canonical Prepared."""

    physical_base: LargeBandPhysicalBaseBinding
    canonical: PrimaryGeometryCanonicalPreparedWitness

    def __post_init__(self) -> None:
        if not isinstance(self.physical_base, LargeBandPhysicalBaseBinding):
            raise TypeError("physical_base must be a LargeBandPhysicalBaseBinding")
        if not isinstance(self.canonical, PrimaryGeometryCanonicalPreparedWitness):
            raise TypeError("canonical must be a PrimaryGeometryCanonicalPreparedWitness")
        if self.canonical.source_key != self.physical_base.source_key:
            raise ValueError("canonical source revision must match the physical-base source")
        if self.canonical.prepared_key != self.physical_base.prepared_key:
            raise ValueError("canonical Prepared identity must match the physical-base Prepared")
        failed = [name for name, ok in self.binding_checks().items() if not ok]
        if failed:
            raise ValueError("uncertified canonical physical-base binding: " + ", ".join(failed))

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.physical_base.box_key

    @property
    def source_key(self) -> tuple[str, str]:
        return self.physical_base.source_key

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.physical_base.prepared_key

    @property
    def sign_pair(self):
        return self.physical_base.sign_pair

    @property
    def canonical_parameter_key(self) -> tuple[int, float, int]:
        return self.canonical.canonical_parameter_key

    def binding_checks(self) -> dict[str, bool]:
        c = self.canonical
        return {
            "physical_base_binding_intact": self.physical_base.physical_base_identity_theorem_certified,
            "source_revision_identity": c.source_key == self.physical_base.source_key,
            "prepared_instance_identity": c.prepared_key == self.physical_base.prepared_key,
            "full_true_cone_witness_certified": c.full_true_cone_witness_certified is True,
            "canonical_prepared_application_certified": c.canonical_prepared_application_certified is True,
            "canonical_phases_application_certified": c.canonical_phases_application_certified is True,
            "canonical_upper_twice_active_right_certified": c.canonical_upper_twice_active_right_certified is True,
            "frequency_uses_canonical_scales_certified": c.frequency_uses_canonical_scales_certified is True,
            "axial_uses_canonical_scales_certified": c.axial_uses_canonical_scales_certified is True,
            "both_signs_share_canonical_phases_certified": c.both_signs_share_canonical_phases_certified is True,
            "theorem_evidence_not_sampled": c.evidence_kind in _ACCEPTED_EVIDENCE,
            "provenance_present": bool(c.provenance),
            "formal_source_pin_intact": (
                c.lean_repository == PINNED_LEAN_REPOSITORY
                and c.lean_commit == PINNED_LEAN_COMMIT
                and c.lean_file == PINNED_LEAN_FILE
                and c.canonical_prepared_def == PINNED_CANONICAL_PREPARED_DEF
                and c.canonical_phases_def == PINNED_CANONICAL_PHASES_DEF
                and c.active_right_def == PINNED_ACTIVE_RIGHT_DEF
                and c.frequency_def == PINNED_FREQUENCY_DEF
                and c.axial_def == PINNED_AXIAL_DEF
                and c.final_slow_base_scales_def == PINNED_FINAL_SLOW_BASE_SCALES_DEF
                and c.full_true_cone_decl == PINNED_FULL_TRUE_CONE
            ),
        }

    @property
    def canonical_physical_base_identity_theorem_certified(self) -> bool:
        return all(self.binding_checks().values())

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_canonical_prepared_application_machine_verified(self) -> bool:
        return False

    @property
    def actual_physical_base_values_materialized(self) -> bool:
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
