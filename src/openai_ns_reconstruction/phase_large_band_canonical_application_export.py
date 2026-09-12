"""Content-addressed handoff for an external canonical Prepared application export.

The canonical scope gate proves only coherence of theorem-facing metadata inside
this repository.  It deliberately does not evaluate Lean's noncomputable
``PrimaryGeometryAssembly.canonicalPrepared``.  The next useful interface is
therefore a strict handoff format for a future Lean/formal runner that exports
one concrete canonical application and the finite active slice attached to it.

This module does not create such an export.  It only validates a supplied
record against an already admitted
``LargeBandCanonicalPhysicalBaseScopeCoverage``.  The record is content
addressed, pinned to the exact formal source, and must reproduce the exact
unsigned boxes, signed labels, Prepared identity, and canonical ``(B,r0,N0)``
tuple.  ``prepared_N`` is recorded as an opaque output of the external formal
runner; Python does not infer or certify its mathematical correctness.

Successful admission remains ``formal-structure``.  A matching SHA-256 proves
payload integrity, not theorem authenticity, so no truth flag for actual Lean
replay, CellIndex enumeration, field materialization, Eqs. (7.9)--(7.11), or
paper-exact velocity is promoted here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Integral
import re

from .phase_large_band_canonical_physical_base import (
    PINNED_CANONICAL_PHASES_DEF,
    PINNED_CANONICAL_PREPARED_DEF,
    PINNED_FULL_TRUE_CONE,
)
from .phase_large_band_canonical_scope_coverage import (
    LargeBandCanonicalPhysicalBaseScopeCoverage,
)
from .phase_large_band_prepared_source import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
)


EXPORT_SCHEMA = "openai-ns/canonical-prepared-application-export/v1"
EXPORT_PRODUCER_KIND = "lean-formal-export"
REQUIRED_APPLICATION_SYMBOLS = (
    PINNED_CANONICAL_PREPARED_DEF,
    PINNED_CANONICAL_PHASES_DEF,
    PINNED_FULL_TRUE_CONE,
)
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _sha256_text(value: object, name: str) -> str:
    text = _nonempty_text(value, name)
    if _SHA256_RE.fullmatch(text) is None:
        raise ValueError(f"{name} must be a lowercase 64-hex SHA-256 digest")
    return text


def _canonical_positive_float_hex(value: object, name: str) -> str:
    text = _nonempty_text(value, name)
    try:
        parsed = float.fromhex(text)
    except ValueError as exc:
        raise ValueError(f"{name} must be a canonical finite positive float.hex() string") from exc
    if not math.isfinite(parsed) or parsed <= 0.0 or parsed.hex() != text:
        raise ValueError(f"{name} must be a canonical finite positive float.hex() string")
    return text


def _box_record(value: object) -> tuple[int, tuple[int, int, int]]:
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError("every unsigned box record must be (ell, (a1,a2,a3))")
    ell = _natural(value[0], "box ell")
    a = value[1]
    if ell == 0:
        raise ValueError("box ell must be positive")
    if not isinstance(a, tuple) or len(a) != 3:
        raise ValueError("every unsigned box record must be (ell, (a1,a2,a3))")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in a):
        raise ValueError("slow-box coordinates must be integers")
    return ell, tuple(int(x) for x in a)


def _label_record(value: object) -> tuple[int, tuple[int, int, int], int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError("every signed label record must be (ell, (a1,a2,a3), sigma)")
    ell, a = _box_record((value[0], value[1]))
    sigma = value[2]
    if isinstance(sigma, bool) or not isinstance(sigma, Integral) or int(sigma) not in (-1, 1):
        raise ValueError("signed label sigma must be -1 or +1")
    return ell, a, int(sigma)


def _digest(payload: dict[str, object]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _scope_box_records(
    coverage: LargeBandCanonicalPhysicalBaseScopeCoverage,
) -> tuple[tuple[int, tuple[int, int, int]], ...]:
    return tuple(sorted(coverage.covered_box_keys, key=lambda item: (item[0], item[1])))


def _scope_label_records(
    coverage: LargeBandCanonicalPhysicalBaseScopeCoverage,
) -> tuple[tuple[int, tuple[int, int, int], int], ...]:
    return tuple(
        sorted(
            ((label.ell, label.a, label.sigma) for label in coverage.covered_signed_labels),
            key=lambda item: (item[0], item[1], item[2]),
        )
    )


def canonical_scope_payload(
    coverage: LargeBandCanonicalPhysicalBaseScopeCoverage,
) -> dict[str, object]:
    """Return the deterministic scope payload expected from a formal export."""

    if not isinstance(coverage, LargeBandCanonicalPhysicalBaseScopeCoverage):
        raise TypeError("coverage must be a LargeBandCanonicalPhysicalBaseScopeCoverage")
    if not coverage.canonical_scope_machine_checked:
        raise ValueError("canonical scope is no longer internally coherent")

    source_id, source_revision, prepared_instance_id = coverage.shared_prepared_key
    B, r0, N0 = coverage.canonical_parameter_key
    boxes = _scope_box_records(coverage)
    labels = _scope_label_records(coverage)
    return {
        "schema": EXPORT_SCHEMA,
        "formal_source": {
            "repository": PINNED_LEAN_REPOSITORY,
            "commit": PINNED_LEAN_COMMIT,
        },
        "source_id": source_id,
        "source_revision": source_revision,
        "prepared_instance_id": prepared_instance_id,
        "canonical_parameters": {
            "B": B,
            "r0_hex": float(r0).hex(),
            "N0": N0,
        },
        "unsigned_boxes": [[ell, list(a)] for ell, a in boxes],
        "signed_labels": [[ell, list(a), sigma] for ell, a, sigma in labels],
        "application_symbols": list(REQUIRED_APPLICATION_SYMBOLS),
    }


def canonical_scope_sha256(
    coverage: LargeBandCanonicalPhysicalBaseScopeCoverage,
) -> str:
    """SHA-256 of the exact canonical scope/identity contract."""

    return _digest(canonical_scope_payload(coverage))


@dataclass(frozen=True)
class CanonicalPreparedApplicationExport:
    """Opaque record supplied by an external formal-export pipeline.

    The class validates shape, canonical encodings and pinned-source identity.
    Exact agreement with a repository scope and the payload digest is checked by
    :class:`LargeBandCanonicalApplicationExportAdmission`.
    """

    schema: str
    producer_kind: str
    formal_repository: str
    formal_commit: str
    source_id: str
    source_revision: str
    prepared_instance_id: str
    prepared_N: int
    B: int
    r0_hex: str
    N0: int
    unsigned_boxes: tuple[tuple[int, tuple[int, int, int]], ...]
    signed_labels: tuple[tuple[int, tuple[int, int, int], int], ...]
    scope_sha256: str
    application_symbols: tuple[str, ...]
    provenance: str
    payload_sha256: str

    def __post_init__(self) -> None:
        if self.schema != EXPORT_SCHEMA:
            raise ValueError("schema must match the canonical export contract")
        if self.producer_kind != EXPORT_PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; sampled/fitted/numeric exports are rejected"
            )
        if self.formal_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("formal_repository must match the pinned Lean source")
        if self.formal_commit != PINNED_LEAN_COMMIT:
            raise ValueError("formal_commit must match the pinned Lean source")

        for name in ("source_id", "source_revision", "prepared_instance_id", "provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))
        object.__setattr__(self, "prepared_N", _natural(self.prepared_N, "prepared_N"))
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "r0_hex", _canonical_positive_float_hex(self.r0_hex, "r0_hex"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "scope_sha256", _sha256_text(self.scope_sha256, "scope_sha256"))
        object.__setattr__(self, "payload_sha256", _sha256_text(self.payload_sha256, "payload_sha256"))

        if not isinstance(self.unsigned_boxes, tuple) or not self.unsigned_boxes:
            raise ValueError("unsigned_boxes must be a nonempty tuple")
        boxes = tuple(_box_record(record) for record in self.unsigned_boxes)
        if len(set(boxes)) != len(boxes):
            raise ValueError("unsigned_boxes must not contain duplicates")
        object.__setattr__(self, "unsigned_boxes", tuple(sorted(boxes)))

        if not isinstance(self.signed_labels, tuple) or not self.signed_labels:
            raise ValueError("signed_labels must be a nonempty tuple")
        labels = tuple(_label_record(record) for record in self.signed_labels)
        if len(set(labels)) != len(labels):
            raise ValueError("signed_labels must not contain duplicates")
        object.__setattr__(self, "signed_labels", tuple(sorted(labels)))

        if not isinstance(self.application_symbols, tuple) or not self.application_symbols:
            raise ValueError("application_symbols must be a nonempty tuple")
        symbols = tuple(_nonempty_text(symbol, "application symbol") for symbol in self.application_symbols)
        if len(set(symbols)) != len(symbols):
            raise ValueError("application_symbols must not contain duplicates")
        object.__setattr__(self, "application_symbols", symbols)


def canonical_application_export_payload(
    export: CanonicalPreparedApplicationExport,
) -> dict[str, object]:
    """Canonical JSON-safe payload, excluding its own payload digest."""

    if not isinstance(export, CanonicalPreparedApplicationExport):
        raise TypeError("export must be a CanonicalPreparedApplicationExport")
    return {
        "schema": export.schema,
        "producer_kind": export.producer_kind,
        "formal_repository": export.formal_repository,
        "formal_commit": export.formal_commit,
        "source_id": export.source_id,
        "source_revision": export.source_revision,
        "prepared_instance_id": export.prepared_instance_id,
        "prepared_N": export.prepared_N,
        "B": export.B,
        "r0_hex": export.r0_hex,
        "N0": export.N0,
        "unsigned_boxes": [[ell, list(a)] for ell, a in export.unsigned_boxes],
        "signed_labels": [
            [ell, list(a), sigma] for ell, a, sigma in export.signed_labels
        ],
        "scope_sha256": export.scope_sha256,
        "application_symbols": list(export.application_symbols),
        "provenance": export.provenance,
    }


def canonical_application_export_sha256(
    export: CanonicalPreparedApplicationExport,
) -> str:
    """Recompute the content digest carried by one export record."""

    return _digest(canonical_application_export_payload(export))


@dataclass(frozen=True)
class LargeBandCanonicalApplicationExportAdmission:
    """Bind one content-addressed external export to one exact canonical scope."""

    canonical_scope: LargeBandCanonicalPhysicalBaseScopeCoverage
    export: CanonicalPreparedApplicationExport

    def __post_init__(self) -> None:
        if not isinstance(self.canonical_scope, LargeBandCanonicalPhysicalBaseScopeCoverage):
            raise TypeError(
                "canonical_scope must be a LargeBandCanonicalPhysicalBaseScopeCoverage"
            )
        if not isinstance(self.export, CanonicalPreparedApplicationExport):
            raise TypeError("export must be a CanonicalPreparedApplicationExport")
        if not self.canonical_scope.canonical_scope_machine_checked:
            raise ValueError("canonical scope is no longer internally coherent")

        expected_scope_digest = canonical_scope_sha256(self.canonical_scope)
        if self.export.scope_sha256 != expected_scope_digest:
            raise ValueError("export scope digest does not match the admitted canonical scope")
        if self.export.payload_sha256 != canonical_application_export_sha256(self.export):
            raise ValueError("export payload SHA-256 does not match its content")

        source_id, source_revision, prepared_instance_id = self.canonical_scope.shared_prepared_key
        if (self.export.source_id, self.export.source_revision, self.export.prepared_instance_id) != (
            source_id,
            source_revision,
            prepared_instance_id,
        ):
            raise ValueError("export Prepared identity does not match the canonical scope")

        B, r0, N0 = self.canonical_scope.canonical_parameter_key
        if (self.export.B, self.export.r0_hex, self.export.N0) != (
            B,
            float(r0).hex(),
            N0,
        ):
            raise ValueError("export canonical (B,r0,N0) parameters do not match the scope")
        if self.export.unsigned_boxes != _scope_box_records(self.canonical_scope):
            raise ValueError("export unsigned box set does not match the canonical scope")
        if self.export.signed_labels != _scope_label_records(self.canonical_scope):
            raise ValueError("export signed label set does not match the canonical scope")
        if self.export.application_symbols != REQUIRED_APPLICATION_SYMBOLS:
            raise ValueError("export application symbols do not match the pinned contract")

    @property
    def export_integrity_machine_checked(self) -> bool:
        return (
            self.export.scope_sha256 == canonical_scope_sha256(self.canonical_scope)
            and self.export.payload_sha256 == canonical_application_export_sha256(self.export)
        )

    @property
    def prepared_N(self) -> int:
        return self.export.prepared_N

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_canonical_prepared_application_machine_verified(self) -> bool:
        return False

    @property
    def actual_cell_index_enumeration_machine_verified(self) -> bool:
        return False

    @property
    def actual_physical_base_values_materialized(self) -> bool:
        return False

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def theorem_applications_machine_replayed(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
