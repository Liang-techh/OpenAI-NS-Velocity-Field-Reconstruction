"""Exact-set coverage of physical-base identities over a declared Section 7 scope.

``phase_large_band_primary_geometry_scope`` ties a theorem-provenanced unsigned
fixed-band manifest to the exact signed family scope, while
``phase_large_band_physical_base`` ties one sign-free slow box to the pinned
``FinalSlowBase.velocity`` identity chain.  Neither layer alone prevents a
caller from supplying physical-base identities for only a strict subset of the
declared active boxes and then describing that subset as family-wide coverage.

This module closes only that bookkeeping seam.  It requires one
``LargeBandPhysicalBaseBinding`` for every unsigned box in an already admitted
``LargeBandPrimaryGeometryScopeAdmission``.  Box keys must match exactly, each
binding must reuse the same base-source revision, duplicates/extras are rejected,
and the union of the two sign labels exposed by every binding must equal the
scope's exact signed label set.

No Lean object is evaluated and no field value is materialized here.  In
particular, the upstream ``CellIndex`` manifest is still an externally
provenanced declaration rather than a Python enumeration of the noncomputable
formal index.  Successful coverage therefore remains ``formal-structure`` and
does not promote Eqs. (7.9)--(7.11) or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass

from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_physical_base import LargeBandPhysicalBaseBinding
from .phase_large_band_primary_geometry_scope import (
    LargeBandPrimaryGeometryScopeAdmission,
)


@dataclass(frozen=True)
class LargeBandPhysicalBaseScopeCoverage:
    """Require exact physical-base identity coverage of one admitted band slice."""

    scope_admission: LargeBandPrimaryGeometryScopeAdmission
    physical_base_bindings: tuple[LargeBandPhysicalBaseBinding, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.scope_admission, LargeBandPrimaryGeometryScopeAdmission):
            raise TypeError(
                "scope_admission must be a LargeBandPrimaryGeometryScopeAdmission"
            )
        if not isinstance(self.physical_base_bindings, tuple):
            raise TypeError("physical_base_bindings must be a tuple")
        if not self.physical_base_bindings:
            raise ValueError("physical_base_bindings must be nonempty")
        if any(
            not isinstance(binding, LargeBandPhysicalBaseBinding)
            for binding in self.physical_base_bindings
        ):
            raise TypeError(
                "every physical-base binding must be a LargeBandPhysicalBaseBinding"
            )

        expected_boxes = self.scope_admission.witness.expected_box_keys
        binding_boxes = [binding.box_key for binding in self.physical_base_bindings]
        if len(set(binding_boxes)) != len(binding_boxes):
            raise ValueError("physical-base coverage must contain exactly one binding per box")
        actual_boxes = frozenset(binding_boxes)
        if actual_boxes != expected_boxes:
            missing = expected_boxes - actual_boxes
            extra = actual_boxes - expected_boxes
            raise ValueError(
                "physical-base box coverage must equal the admitted primary-geometry scope; "
                f"missing={len(missing)}, extra={len(extra)}"
            )

        expected_source = self.scope_admission.scope.source_key
        for binding in self.physical_base_bindings:
            if binding.source_key != expected_source:
                raise ValueError(
                    "every physical-base binding must reuse the admitted scope source revision"
                )
            if not binding.physical_base_identity_theorem_certified:
                raise ValueError("physical-base binding is no longer internally certified")

        if self.covered_signed_labels != self.scope_admission.scope.expected_labels:
            raise ValueError(
                "physical-base sign-pair coverage must equal the admitted signed family scope"
            )

    @property
    def covered_box_keys(self) -> frozenset[tuple[int, tuple[int, int, int]]]:
        return frozenset(binding.box_key for binding in self.physical_base_bindings)

    @property
    def covered_signed_labels(self) -> frozenset[AsymptoticSlowLabel]:
        return frozenset(
            label
            for binding in self.physical_base_bindings
            for label in binding.sign_pair
        )

    @property
    def physical_binding_count(self) -> int:
        return len(self.physical_base_bindings)

    @property
    def declared_scope_physical_base_coverage_machine_checked(self) -> bool:
        """Whether Python checked exact-set composition relative to the supplied scope."""

        return (
            self.covered_box_keys == self.scope_admission.witness.expected_box_keys
            and self.covered_signed_labels == self.scope_admission.scope.expected_labels
        )

    @property
    def status(self) -> str:
        return "formal-structure"

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
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
