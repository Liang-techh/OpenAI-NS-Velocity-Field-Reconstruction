"""Fail-closed bridge from a formal ``Prepared`` source to the physical slow base.

The previous large-band adapter identifies a sign-free slow box with the pinned
``PrimaryGeometryAssembly.Prepared`` construction, but it deliberately stops at
constructor identity.  The pinned formal source contains a stronger field-level
chain on every active carrier:

``construction_frequency`` / ``construction_axial`` identify the two phase
base fields with ``frequency`` / ``axial`` at ``cellBand L``; and
``frequency_eq_physical`` / ``axial_eq_physical`` identify those fields with
scaled components 1 and 2 of the same ``FinalSlowBase.velocity``.  Carrier-time
positivity and ``Prepared.radius_pos`` discharge the positivity hypotheses on
the active carrier.

This module records an externally theorem-certified application of exactly that
chain for one admitted large-band box.  It evaluates no Lean object, samples no
field, and constructs no surrogate background.  Therefore successful admission
remains ``formal-structure`` until a machine-linked theorem/export supplies the
actual applications and values.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .phase_large_band_prepared_source import (
    LargeBandPreparedSourceBinding,
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_FILE,
    PINNED_LEAN_REPOSITORY,
)


PINNED_CONSTRUCTION_FREQUENCY = "PrimaryGeometryAssembly.construction_frequency"
PINNED_CONSTRUCTION_AXIAL = "PrimaryGeometryAssembly.construction_axial"
PINNED_CARRIER_TIME_POSITIVE = "PrimaryGeometryAssembly.carrier_time_positive"
PINNED_PREPARED_RADIUS_POS = "PrimaryGeometryAssembly.Prepared.radius_pos"
PINNED_FREQUENCY_EQ_PHYSICAL = "PrimaryGeometryAssembly.frequency_eq_physical"
PINNED_AXIAL_EQ_PHYSICAL = "PrimaryGeometryAssembly.axial_eq_physical"
PINNED_FINAL_SLOW_BASE_VELOCITY = "FinalSlowBase.velocity"
PINNED_CELL_BAND = "BaseChartJets.cellBand"

_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _box3(value: object) -> tuple[int, int, int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError("a must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in value):
        raise ValueError("a must be a length-three integer tuple")
    return tuple(int(x) for x in value)


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
class PrimaryGeometryPhysicalBaseWitness:
    """Theorem-facing identity of one active box with ``FinalSlowBase.velocity``.

    The witness is intentionally application metadata.  The certification flags
    must come from a theorem/export layer; Python only checks that all pieces of
    the pinned composition refer to one box, one source revision, and one
    selected ``Prepared`` object.
    """

    ell: int
    a: tuple[int, int, int]
    source_id: str
    source_revision: str
    prepared_instance_id: str
    evidence_kind: str
    provenance: str
    local_base_same_frequency_axial_certified: bool
    construction_frequency_certified: bool
    construction_axial_certified: bool
    carrier_time_positive_certified: bool
    prepared_radius_positive_certified: bool
    frequency_eq_physical_certified: bool
    axial_eq_physical_certified: bool
    both_signs_share_physical_base_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_LEAN_FILE
    construction_frequency_theorem: str = PINNED_CONSTRUCTION_FREQUENCY
    construction_axial_theorem: str = PINNED_CONSTRUCTION_AXIAL
    carrier_time_positive_theorem: str = PINNED_CARRIER_TIME_POSITIVE
    prepared_radius_pos_field: str = PINNED_PREPARED_RADIUS_POS
    frequency_eq_physical_theorem: str = PINNED_FREQUENCY_EQ_PHYSICAL
    axial_eq_physical_theorem: str = PINNED_AXIAL_EQ_PHYSICAL
    physical_velocity_definition: str = PINNED_FINAL_SLOW_BASE_VELOCITY
    cell_band_definition: str = PINNED_CELL_BAND

    def __post_init__(self) -> None:
        object.__setattr__(self, "ell", _positive_integer(self.ell, "ell"))
        object.__setattr__(self, "a", _box3(self.a))
        for name in ("source_id", "source_revision", "prepared_instance_id", "provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))
        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        for name in (
            "local_base_same_frequency_axial_certified",
            "construction_frequency_certified",
            "construction_axial_certified",
            "carrier_time_positive_certified",
            "prepared_radius_positive_certified",
            "frequency_eq_physical_certified",
            "axial_eq_physical_certified",
            "both_signs_share_physical_base_certified",
        ):
            _strict_true(getattr(self, name), name)

        pins = (
            ("lean_repository", PINNED_LEAN_REPOSITORY),
            ("lean_commit", PINNED_LEAN_COMMIT),
            ("lean_file", PINNED_LEAN_FILE),
            ("construction_frequency_theorem", PINNED_CONSTRUCTION_FREQUENCY),
            ("construction_axial_theorem", PINNED_CONSTRUCTION_AXIAL),
            ("carrier_time_positive_theorem", PINNED_CARRIER_TIME_POSITIVE),
            ("prepared_radius_pos_field", PINNED_PREPARED_RADIUS_POS),
            ("frequency_eq_physical_theorem", PINNED_FREQUENCY_EQ_PHYSICAL),
            ("axial_eq_physical_theorem", PINNED_AXIAL_EQ_PHYSICAL),
            ("physical_velocity_definition", PINNED_FINAL_SLOW_BASE_VELOCITY),
            ("cell_band_definition", PINNED_CELL_BAND),
        )
        for name, expected in pins:
            object.__setattr__(self, name, _pinned(getattr(self, name), expected, name))

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.ell, self.a

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.source_id, self.source_revision, self.prepared_instance_id


@dataclass(frozen=True)
class LargeBandPhysicalBaseBinding:
    """Bind one admitted slow box to the pinned physical-base identity chain."""

    prepared_source: LargeBandPreparedSourceBinding
    witness: PrimaryGeometryPhysicalBaseWitness

    def __post_init__(self) -> None:
        if not isinstance(self.prepared_source, LargeBandPreparedSourceBinding):
            raise TypeError("prepared_source must be a LargeBandPreparedSourceBinding")
        if not isinstance(self.witness, PrimaryGeometryPhysicalBaseWitness):
            raise TypeError("witness must be a PrimaryGeometryPhysicalBaseWitness")
        if self.witness.box_key != self.prepared_source.base_source.box_key:
            raise ValueError("physical-base box key must match the admitted prepared slow box")
        if self.witness.source_key != self.prepared_source.source_key:
            raise ValueError("physical-base source revision must match the admitted prepared source")
        if self.witness.prepared_key != self.prepared_source.prepared_key:
            raise ValueError("physical-base Prepared identity must match the selected prepared source")
        failed = [name for name, ok in self.binding_checks().items() if not ok]
        if failed:
            raise ValueError("uncertified physical-base binding: " + ", ".join(failed))

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.witness.box_key

    @property
    def source_key(self) -> tuple[str, str]:
        return self.witness.source_key

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.witness.prepared_key

    @property
    def sign_pair(self):
        return self.prepared_source.sign_pair

    def binding_checks(self) -> dict[str, bool]:
        return {
            "prepared_source_binding_intact": all(self.prepared_source.binding_checks().values()),
            "box_identity": self.witness.box_key == self.prepared_source.base_source.box_key,
            "source_revision_identity": self.witness.source_key == self.prepared_source.source_key,
            "prepared_instance_identity": self.witness.prepared_key == self.prepared_source.prepared_key,
            "local_base_same_frequency_axial_certified": self.witness.local_base_same_frequency_axial_certified is True,
            "construction_frequency_certified": self.witness.construction_frequency_certified is True,
            "construction_axial_certified": self.witness.construction_axial_certified is True,
            "carrier_time_positive_certified": self.witness.carrier_time_positive_certified is True,
            "prepared_radius_positive_certified": self.witness.prepared_radius_positive_certified is True,
            "frequency_eq_physical_certified": self.witness.frequency_eq_physical_certified is True,
            "axial_eq_physical_certified": self.witness.axial_eq_physical_certified is True,
            "both_signs_share_physical_base_certified": self.witness.both_signs_share_physical_base_certified is True,
            "theorem_evidence_not_sampled": self.witness.evidence_kind in _ACCEPTED_EVIDENCE,
            "provenance_present": bool(self.witness.provenance),
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def physical_base_identity_theorem_certified(self) -> bool:
        return all(self.binding_checks().values())

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
