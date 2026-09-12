"""Bind a declared Section 7 family scope to pinned primary-geometry indexing.

The finite-scope coverage gate in :mod:`phase_large_band_family_coverage`
intentionally treats its active signed label tuple as externally supplied.
That prevents sparse theorem witnesses from masquerading as uniform coverage,
but it does not identify the tuple with the formal primary-geometry index used
by the paper reconstruction.

At the pinned Lean revision, ``BaseChartJets.CellIndex`` is the subtype of
actual positive active labels above the fixed lower band, and
``PrimaryGeometryAssembly.cellDomain`` uses that index for the carrier domain.
``PrimaryGeometryAssembly.family`` then builds the two ``Fin 2`` phase-family
branches on the same index.  This module provides a fail-closed interface for
an *externally proved/exported* fixed-band slice of those definitions.

Python does not enumerate ``CellIndex`` and does not infer a manuscript active
set.  A theorem-level witness must supply the complete unsigned box manifest
for one band and certify its projection from the pinned formal index.  The
admission layer only performs the deterministic part: it duplicates every
unsigned box into both signs and requires exact equality with an existing
``LargeBandActiveFamilyScopeWitness``.  Samples, fitted manifests, and numeric
scans are rejected.

Consequently a successful admission is still ``formal-structure``.  It is not
a proof that the supplied manifest is the paper's actual active family, does
not materialize Proposition 5.5 base fields, and does not promote Eqs.
(7.9)--(7.11) or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .phase_large_band_family_coverage import LargeBandActiveFamilyScopeWitness
from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_phase_estimates import PINNED_LEAN_COMMIT, PINNED_LEAN_REPOSITORY


PINNED_CELL_INDEX_DEFINITION = "BaseChartJets.CellIndex"
PINNED_CELL_DOMAIN_DEFINITION = "PrimaryGeometryAssembly.cellDomain"
PINNED_FAMILY_DEFINITION = "PrimaryGeometryAssembly.family"

_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _box3(value: object) -> tuple[int, int, int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError("every unsigned box must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in value):
        raise ValueError("every unsigned box must be a length-three integer tuple")
    return tuple(int(x) for x in value)


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class PrimaryGeometryBandSliceWitness:
    """Theorem-facing export of one fixed-band unsigned ``CellIndex`` slice.

    ``unsigned_boxes`` is data, not proof.  Its completeness and its projection
    from the pinned formal definitions are represented by the required theorem
    certifications below.  The Python verifier then independently reconstructs
    the signed label set from this manifest and compares it with the declared
    Section 7 scope.

    No convention ``Fin 2 = {0,1} -> sigma = {-1,+1}`` is hard-coded.  Only the
    theorem-certified fact that the two formal branches correspond bijectively
    to the two Section 6 signs is consumed.
    """

    scope_id: str
    ell: int
    source_id: str
    source_revision: str
    unsigned_boxes: tuple[tuple[int, int, int], ...]
    lean_repository: str
    lean_commit: str
    cell_index_definition: str
    cell_domain_definition: str
    family_definition: str
    evidence_kind: str
    provenance: str
    complete_unsigned_box_manifest_certified: bool
    cell_index_band_slice_certified: bool
    band_and_cell_projection_certified: bool
    cell_domain_carrier_identity_certified: bool
    family_index_identity_certified: bool
    fin2_sign_duplication_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "scope_id", _nonempty_text(self.scope_id, "scope_id"))
        object.__setattr__(self, "ell", _positive_integer(self.ell, "ell"))
        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self, "source_revision", _nonempty_text(self.source_revision, "source_revision")
        )
        object.__setattr__(self, "provenance", _nonempty_text(self.provenance, "provenance"))

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted/numeric-scan manifests are rejected"
            )
        if self.lean_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("lean_repository must match the pinned formal source")
        if self.lean_commit != PINNED_LEAN_COMMIT:
            raise ValueError("lean_commit must match the pinned formal source")
        if self.cell_index_definition != PINNED_CELL_INDEX_DEFINITION:
            raise ValueError("cell_index_definition must match BaseChartJets.CellIndex")
        if self.cell_domain_definition != PINNED_CELL_DOMAIN_DEFINITION:
            raise ValueError("cell_domain_definition must match PrimaryGeometryAssembly.cellDomain")
        if self.family_definition != PINNED_FAMILY_DEFINITION:
            raise ValueError("family_definition must match PrimaryGeometryAssembly.family")

        if not isinstance(self.unsigned_boxes, tuple) or not self.unsigned_boxes:
            raise ValueError("unsigned_boxes must be a nonempty tuple")
        boxes = tuple(_box3(box) for box in self.unsigned_boxes)
        if len(set(boxes)) != len(boxes):
            raise ValueError("unsigned_boxes must be unique")
        object.__setattr__(self, "unsigned_boxes", boxes)

        for name in (
            "complete_unsigned_box_manifest_certified",
            "cell_index_band_slice_certified",
            "band_and_cell_projection_certified",
            "cell_domain_carrier_identity_certified",
            "family_index_identity_certified",
            "fin2_sign_duplication_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def expected_box_keys(self) -> frozenset[tuple[int, tuple[int, int, int]]]:
        return frozenset((self.ell, box) for box in self.unsigned_boxes)

    @property
    def expected_signed_labels(self) -> frozenset[AsymptoticSlowLabel]:
        return frozenset(
            AsymptoticSlowLabel(self.ell, box, sigma)
            for box in self.unsigned_boxes
            for sigma in (-1, 1)
        )


@dataclass(frozen=True)
class LargeBandPrimaryGeometryScopeAdmission:
    """Exact-set bridge from a pinned formal band-slice manifest to #164 scope."""

    scope: LargeBandActiveFamilyScopeWitness
    witness: PrimaryGeometryBandSliceWitness

    def __post_init__(self) -> None:
        if not isinstance(self.scope, LargeBandActiveFamilyScopeWitness):
            raise TypeError("scope must be a LargeBandActiveFamilyScopeWitness")
        if not isinstance(self.witness, PrimaryGeometryBandSliceWitness):
            raise TypeError("witness must be a PrimaryGeometryBandSliceWitness")

        if self.scope.scope_id != self.witness.scope_id:
            raise ValueError("scope_id must match the primary-geometry band-slice witness")
        if self.scope.ell != self.witness.ell:
            raise ValueError("ell must match the primary-geometry band-slice witness")
        if self.scope.source_key != self.witness.source_key:
            raise ValueError("base source revision must match the primary-geometry band-slice witness")
        if self.scope.box_keys != self.witness.expected_box_keys:
            raise ValueError("unsigned box set must equal the declared active scope box set")
        if self.scope.expected_labels != self.witness.expected_signed_labels:
            raise ValueError("signed label set must equal the theorem-provenanced primary-geometry band slice")

    @property
    def unsigned_box_count(self) -> int:
        return len(self.witness.unsigned_boxes)

    @property
    def signed_family_count(self) -> int:
        return len(self.scope.labels)

    @property
    def manifest_matches_scope(self) -> bool:
        return (
            self.scope.box_keys == self.witness.expected_box_keys
            and self.scope.expected_labels == self.witness.expected_signed_labels
        )

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_cell_index_enumeration_machine_verified(self) -> bool:
        return False

    @property
    def actual_active_family_scope_machine_verified(self) -> bool:
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
