"""Exact-scope coverage for the pinned canonical Prepared object.

``phase_large_band_canonical_physical_base`` certifies, one slow box at a
time, that an admitted ``Prepared -> FinalSlowBase.velocity`` binding uses the
pinned ``PrimaryGeometryAssembly.canonicalPrepared`` / ``canonicalPhases``
path.  The formal constructor is stronger than a collection of unrelated
per-box applications: ``canonicalPrepared`` returns one ``Prepared`` object,
and that object's ``base`` field is indexed over every ``L : Index W N``.

Without an additional scope check, callers could therefore provide a distinct
canonical-looking ``Prepared`` identity for each box and describe the set as a
single canonical primary-geometry family.  This module closes only that seam.
It requires exact box coverage of an existing physical-base scope, requires
each canonical binding to wrap the exact physical binding already admitted for
that box, and requires one shared Prepared identity and one shared ``(B,r0,N0)``
canonical constructor tuple across the whole scope.

No Lean term is evaluated here.  The supplied canonical applications remain
theorem-facing provenance metadata, so successful coverage is still only
``formal-structure`` and does not promote Eqs. (7.9)--(7.11) or paper-exact
velocity.
"""
from __future__ import annotations

from dataclasses import dataclass

from .phase_large_band_canonical_physical_base import (
    LargeBandCanonicalPhysicalBaseBinding,
)
from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_physical_base_coverage import (
    LargeBandPhysicalBaseScopeCoverage,
)


@dataclass(frozen=True)
class LargeBandCanonicalPhysicalBaseScopeCoverage:
    """Require one globally shared canonical Prepared across an exact band scope."""

    physical_base_coverage: LargeBandPhysicalBaseScopeCoverage
    canonical_bindings: tuple[LargeBandCanonicalPhysicalBaseBinding, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.physical_base_coverage, LargeBandPhysicalBaseScopeCoverage):
            raise TypeError(
                "physical_base_coverage must be a LargeBandPhysicalBaseScopeCoverage"
            )
        if not isinstance(self.canonical_bindings, tuple):
            raise TypeError("canonical_bindings must be a tuple")
        if not self.canonical_bindings:
            raise ValueError("canonical_bindings must be nonempty")
        if any(
            not isinstance(binding, LargeBandCanonicalPhysicalBaseBinding)
            for binding in self.canonical_bindings
        ):
            raise TypeError(
                "every canonical binding must be a LargeBandCanonicalPhysicalBaseBinding"
            )
        if not self.physical_base_coverage.declared_scope_physical_base_coverage_machine_checked:
            raise ValueError("physical-base exact-set coverage is no longer internally complete")

        expected_by_box = {
            binding.box_key: binding
            for binding in self.physical_base_coverage.physical_base_bindings
        }
        canonical_boxes = [binding.box_key for binding in self.canonical_bindings]
        if len(set(canonical_boxes)) != len(canonical_boxes):
            raise ValueError("canonical coverage must contain exactly one binding per box")
        if frozenset(canonical_boxes) != self.physical_base_coverage.covered_box_keys:
            missing = self.physical_base_coverage.covered_box_keys - frozenset(canonical_boxes)
            extra = frozenset(canonical_boxes) - self.physical_base_coverage.covered_box_keys
            raise ValueError(
                "canonical box coverage must equal the admitted physical-base scope; "
                f"missing={len(missing)}, extra={len(extra)}"
            )

        for binding in self.canonical_bindings:
            expected = expected_by_box[binding.box_key]
            if binding.physical_base != expected:
                raise ValueError(
                    "canonical binding must wrap the exact physical-base binding admitted for its box"
                )
            if not binding.canonical_physical_base_identity_theorem_certified:
                raise ValueError("canonical physical-base binding is no longer internally certified")

        prepared_keys = {binding.prepared_key for binding in self.canonical_bindings}
        if len(prepared_keys) != 1:
            raise ValueError(
                "all boxes in one canonical scope must share one canonical Prepared instance"
            )
        parameter_keys = {
            binding.canonical_parameter_key for binding in self.canonical_bindings
        }
        if len(parameter_keys) != 1:
            raise ValueError(
                "all boxes in one canonical scope must share the same canonical (B,r0,N0) tuple"
            )

        if self.covered_signed_labels != self.physical_base_coverage.covered_signed_labels:
            raise ValueError(
                "canonical sign-pair coverage must equal the admitted physical-base signed scope"
            )

    @property
    def covered_box_keys(self) -> frozenset[tuple[int, tuple[int, int, int]]]:
        return frozenset(binding.box_key for binding in self.canonical_bindings)

    @property
    def covered_signed_labels(self) -> frozenset[AsymptoticSlowLabel]:
        return frozenset(
            label for binding in self.canonical_bindings for label in binding.sign_pair
        )

    @property
    def shared_prepared_key(self) -> tuple[str, str, str]:
        return self.canonical_bindings[0].prepared_key

    @property
    def canonical_parameter_key(self) -> tuple[int, float, int]:
        return self.canonical_bindings[0].canonical_parameter_key

    @property
    def canonical_binding_count(self) -> int:
        return len(self.canonical_bindings)

    @property
    def canonical_scope_machine_checked(self) -> bool:
        """Mechanical coherence relative to the supplied theorem-facing scope."""

        return (
            self.covered_box_keys == self.physical_base_coverage.covered_box_keys
            and self.covered_signed_labels
            == self.physical_base_coverage.covered_signed_labels
            and len({binding.prepared_key for binding in self.canonical_bindings}) == 1
            and len(
                {binding.canonical_parameter_key for binding in self.canonical_bindings}
            )
            == 1
        )

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
