"""Bind canonical Section 6 slow-box labels to actual signed Section 7 applications.

The repository has two deliberately separate formal-structure chains:

* ``LargeBandCanonicalApplicationExportAdmission`` content-addresses one
  canonical ``PrimaryGeometryAssembly.canonicalPrepared`` scope and its exact
  finite signed label roster; and
* ``ActualSignedInvariantCommonCurlAdmission`` composes cycle-invariant jets
  with the pinned actual signed common-curl/divergence theorem.

This module closes only the *identity handoff* between those chains.  A future
formal exporter may state that one opaque Lean ``SignedLabel`` is the image of
one canonical ``(ell, a, sigma)`` label, while keeping Lean's ``SignedLabel``
encoding opaque to Python.  The gate requires exact scope coverage, the same
``B/N0`` and selected Prepared ``N``, the same pinned formal revision, and an
already admitted invariant/common-curl application for every mapped label.

No mapping theorem is invented here.  The link witness is explicitly an
external ``lean-formal-export`` record.  Its booleans are fail-closed exporter
facts, not replayed proof terms.  Therefore successful admission remains
``formal-structure`` and does not materialize a wave, prove a pointwise/global
divergence statement beyond the admitted theorem metadata, or make the
paper-exact velocity available.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import re

from .actual_signed_common_curl import (
    PINNED_COMMON_CURL_THEOREM,
    PINNED_COPIES_DEFINITION,
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
)
from .actual_signed_common_curl_invariant import (
    ActualSignedInvariantCommonCurlAdmission,
    PINNED_FULL_REQUEST,
    PINNED_JETS_THEOREM,
    PINNED_PHASE_CELL,
)
from .phase_large_band_canonical_application_export import (
    LargeBandCanonicalApplicationExportAdmission,
)
from .phase_large_band_canonical_physical_base import PINNED_CANONICAL_PHASES_DEF


LINK_SCHEMA = "openai-ns/canonical-to-actual-signed-label-link/v1"
LINK_PRODUCER_KIND = "lean-formal-export"
REQUIRED_LINK_SYMBOLS = (
    PINNED_CANONICAL_PHASES_DEF,
    PINNED_COPIES_DEFINITION,
    PINNED_PHASE_CELL,
    PINNED_FULL_REQUEST,
    PINNED_JETS_THEOREM,
    PINNED_COMMON_CURL_THEOREM,
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


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


def _sha256(value: object, name: str) -> str:
    text = _nonempty_text(value, name)
    if _SHA256_RE.fullmatch(text) is None:
        raise ValueError(f"{name} must be a lowercase 64-hex SHA-256 digest")
    return text


def _slow_box(value: object) -> tuple[int, int, int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError("a must be a three-integer slow-box coordinate")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in value):
        raise ValueError("a must be a three-integer slow-box coordinate")
    return tuple(int(x) for x in value)


@dataclass(frozen=True)
class CanonicalActualSignedLabelLinkWitness:
    """External formal-export identity for one canonical/actual signed label.

    ``signed_label_repr`` remains an opaque stable identifier supplied by the
    formal exporter.  Python deliberately does not encode or reconstruct Lean's
    ``SignedLabel B N0``.
    """

    ell: int
    a: tuple[int, int, int]
    sigma: int
    signed_label_repr: str
    scope_sha256: str
    prepared_N: int
    application_id: str
    producer_kind: str
    provenance: str
    canonical_family_label_identity_certified: bool
    signed_label_export_identity_certified: bool
    canonical_phase_cell_identity_certified: bool
    same_prepared_application_certified: bool
    formal_export_application_certified: bool
    schema: str = LINK_SCHEMA
    formal_repository: str = PINNED_LEAN_REPOSITORY
    formal_commit: str = PINNED_LEAN_COMMIT
    dependency_symbols: tuple[str, ...] = REQUIRED_LINK_SYMBOLS

    def __post_init__(self) -> None:
        ell = _natural(self.ell, "ell")
        if ell == 0:
            raise ValueError("ell must be positive")
        object.__setattr__(self, "ell", ell)
        object.__setattr__(self, "a", _slow_box(self.a))
        if isinstance(self.sigma, bool) or not isinstance(self.sigma, Integral) or int(self.sigma) not in (-1, 1):
            raise ValueError("sigma must be -1 or +1")
        object.__setattr__(self, "sigma", int(self.sigma))
        object.__setattr__(self, "signed_label_repr", _nonempty_text(self.signed_label_repr, "signed_label_repr"))
        object.__setattr__(self, "scope_sha256", _sha256(self.scope_sha256, "scope_sha256"))
        object.__setattr__(self, "prepared_N", _natural(self.prepared_N, "prepared_N"))
        object.__setattr__(self, "application_id", _nonempty_text(self.application_id, "application_id"))
        object.__setattr__(self, "provenance", _nonempty_text(self.provenance, "provenance"))

        if self.schema != LINK_SCHEMA:
            raise ValueError("schema must match the canonical-to-actual signed link contract")
        if self.producer_kind != LINK_PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.formal_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("formal_repository must match the pinned formal source exactly")
        if self.formal_commit != PINNED_LEAN_COMMIT:
            raise ValueError("formal_commit must match the pinned formal source exactly")

        for name in (
            "canonical_family_label_identity_certified",
            "signed_label_export_identity_certified",
            "canonical_phase_cell_identity_certified",
            "same_prepared_application_certified",
            "formal_export_application_certified",
        ):
            _strict_true(getattr(self, name), name)

        if not isinstance(self.dependency_symbols, tuple):
            raise ValueError("dependency_symbols must be the pinned formal-symbol tuple")
        symbols = tuple(_nonempty_text(x, "dependency symbol") for x in self.dependency_symbols)
        if symbols != REQUIRED_LINK_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned formal source exactly")
        object.__setattr__(self, "dependency_symbols", symbols)

    @property
    def canonical_label(self) -> tuple[int, tuple[int, int, int], int]:
        return self.ell, self.a, self.sigma


@dataclass(frozen=True)
class CanonicalSignedInvariantCommonCurlBinding:
    """Bind one formal-export label link to one admitted actual application."""

    link: CanonicalActualSignedLabelLinkWitness
    actual_application: ActualSignedInvariantCommonCurlAdmission

    def __post_init__(self) -> None:
        if not isinstance(self.link, CanonicalActualSignedLabelLinkWitness):
            raise TypeError("link must be a CanonicalActualSignedLabelLinkWitness")
        if not isinstance(self.actual_application, ActualSignedInvariantCommonCurlAdmission):
            raise TypeError("actual_application must be an ActualSignedInvariantCommonCurlAdmission")
        if self.actual_application.common_curl.witness.signed_label_repr != self.link.signed_label_repr:
            raise ValueError("opaque SignedLabel identity does not match the admitted actual application")
        if not self.actual_application.invariant_supplies_common_curl_jets:
            raise ValueError("actual application is no longer an admitted invariant/common-curl composition")

    @property
    def canonical_label(self) -> tuple[int, tuple[int, int, int], int]:
        return self.link.canonical_label


@dataclass(frozen=True)
class LargeBandCanonicalSignedCommonCurlScopeAdmission:
    """Exact-scope identity bridge from canonical labels to actual signed uses."""

    canonical_export: LargeBandCanonicalApplicationExportAdmission
    bindings: tuple[CanonicalSignedInvariantCommonCurlBinding, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.canonical_export, LargeBandCanonicalApplicationExportAdmission):
            raise TypeError("canonical_export must be a LargeBandCanonicalApplicationExportAdmission")
        if not self.canonical_export.export_integrity_machine_checked:
            raise ValueError("canonical export is no longer internally coherent")
        if not isinstance(self.bindings, tuple) or not self.bindings:
            raise ValueError("bindings must be a nonempty tuple")
        if any(not isinstance(item, CanonicalSignedInvariantCommonCurlBinding) for item in self.bindings):
            raise TypeError("every binding must be a CanonicalSignedInvariantCommonCurlBinding")

        export = self.canonical_export.export
        expected_labels = set(export.signed_labels)
        labels = tuple(item.canonical_label for item in self.bindings)
        if len(set(labels)) != len(labels):
            raise ValueError("canonical signed labels must not be duplicated")
        if set(labels) != expected_labels:
            raise ValueError("actual signed bindings must cover the canonical signed label scope exactly")

        opaque_labels = tuple(item.link.signed_label_repr for item in self.bindings)
        if len(set(opaque_labels)) != len(opaque_labels):
            raise ValueError("one opaque Lean SignedLabel cannot represent multiple canonical labels")

        application_ids = tuple(item.link.application_id for item in self.bindings)
        if len(set(application_ids)) != len(application_ids):
            raise ValueError("formal-export application_id values must be unique across the scope")

        for item in self.bindings:
            link = item.link
            actual = item.actual_application
            if link.scope_sha256 != export.scope_sha256:
                raise ValueError("signed-label link scope digest does not match the canonical export")
            if link.prepared_N != export.prepared_N:
                raise ValueError("signed-label link Prepared N does not match the canonical export")
            if (actual.jets.B, actual.jets.N0) != (export.B, export.N0):
                raise ValueError("actual signed application B/N0 does not match the canonical scope")
            if actual.jets.lean_repository != export.formal_repository:
                raise ValueError("actual signed application repository does not match the canonical export")
            if actual.jets.lean_commit != export.formal_commit:
                raise ValueError("actual signed application commit does not match the canonical export")

    @property
    def covered_canonical_labels(self) -> frozenset[tuple[int, tuple[int, int, int], int]]:
        return frozenset(item.canonical_label for item in self.bindings)

    @property
    def canonical_to_actual_signed_scope_machine_checked(self) -> bool:
        return self.covered_canonical_labels == frozenset(self.canonical_export.export.signed_labels)

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_signed_label_export_machine_verified(self) -> bool:
        return False

    @property
    def canonical_prepared_application_machine_replayed(self) -> bool:
        return False

    @property
    def actual_cycle_invariant_machine_verified(self) -> bool:
        return False

    @property
    def common_curl_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_signed_wave_values_materialized(self) -> bool:
        return False

    @property
    def paper_exact_divergence_free_wave_available(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
