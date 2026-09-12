"""Coherent physical-base/estimate coverage for one Section 7 large-band scope.

Two fail-closed coverage layers already exist on the Section 6->7 path:

* :mod:`phase_large_band_physical_base_coverage` requires one pinned
  ``Prepared -> FinalSlowBase.velocity`` identity binding for every unsigned
  box in an admitted primary-geometry band slice; and
* :mod:`phase_large_band_family_coverage` requires the pinned phase,
  coordinate-frame, and damping theorem admissions to cover one exact signed
  family scope.

Those layers can still be instantiated independently.  Without a composition
check, a caller could accidentally pair physical-base coverage for one scope
with theorem-estimate coverage for a different (but superficially similar)
scope and then describe the pair as a uniform Section 7 result.

This module closes only that seam.  It requires both coverage objects to share
exactly the same theorem-facing scope witness, box set, signed-label set and
base-source revision.  For every signed family admission it also checks that
the upstream sign-free base-source binding has the same box, source revision
and ``M`` as the physical ``Prepared`` binding for that box.

No field is evaluated here.  In particular, this does not enumerate the
noncomputable formal ``CellIndex``, materialize ``Prepared`` or
``FinalSlowBase.velocity``, or independently replay the Lean theorem
applications.  Successful composition therefore remains ``formal-structure``
and does not promote Eqs. (7.9)--(7.11) or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass

from .phase_large_band_family_coverage import LargeBandFamilyEstimateCoverage
from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_physical_base import LargeBandPhysicalBaseBinding
from .phase_large_band_physical_base_coverage import (
    LargeBandPhysicalBaseScopeCoverage,
)


@dataclass(frozen=True)
class LargeBandPhysicalEstimateScopeCoverage:
    """Compose physical-base and theorem-estimate coverage on one exact scope."""

    physical_base_coverage: LargeBandPhysicalBaseScopeCoverage
    family_estimate_coverage: LargeBandFamilyEstimateCoverage

    def __post_init__(self) -> None:
        if not isinstance(
            self.physical_base_coverage, LargeBandPhysicalBaseScopeCoverage
        ):
            raise TypeError(
                "physical_base_coverage must be a LargeBandPhysicalBaseScopeCoverage"
            )
        if not isinstance(self.family_estimate_coverage, LargeBandFamilyEstimateCoverage):
            raise TypeError(
                "family_estimate_coverage must be a LargeBandFamilyEstimateCoverage"
            )

        physical_scope = self.physical_base_coverage.scope_admission.scope
        estimate_scope = self.family_estimate_coverage.scope
        if physical_scope != estimate_scope:
            raise ValueError(
                "physical-base and estimate coverage must share the same exact active scope witness"
            )
        if not self.physical_base_coverage.declared_scope_physical_base_coverage_machine_checked:
            raise ValueError("physical-base exact-set coverage is no longer internally complete")
        if not self.family_estimate_coverage.scope_coverage_complete:
            raise ValueError("family-estimate exact-set coverage is no longer internally complete")

        if self.physical_base_coverage.covered_signed_labels != self.family_estimate_coverage.covered_labels:
            raise ValueError(
                "physical-base and theorem-estimate signed-label coverage must be identical"
            )
        if self.physical_base_coverage.covered_box_keys != estimate_scope.box_keys:
            raise ValueError(
                "physical-base and theorem-estimate unsigned box coverage must be identical"
            )

        for phase in self.family_estimate_coverage.phase_admissions:
            binding = self.physical_binding_for_label(phase.label)
            family_source = phase.family.binding
            if family_source.box_key != binding.box_key:
                raise ValueError(
                    "estimate family base-source box must match the physical Prepared box"
                )
            if family_source.source_key != binding.source_key:
                raise ValueError(
                    "estimate family must reuse the physical-base source revision"
                )
            if family_source.source.M != binding.prepared_source.base_source.source.M:
                raise ValueError(
                    "estimate family and physical Prepared binding must use the same LocalBase M"
                )
            if phase.label not in binding.sign_pair:
                raise ValueError(
                    "signed estimate family must belong to the physical binding sign pair"
                )
            if not all(family_source.binding_checks().values()):
                raise ValueError("estimate family base-source binding is no longer certified")
            if not binding.physical_base_identity_theorem_certified:
                raise ValueError("physical-base binding is no longer theorem-certified")

    def physical_binding_for_label(
        self, label: AsymptoticSlowLabel
    ) -> LargeBandPhysicalBaseBinding:
        """Return the unique sign-free physical binding for ``label``'s slow box."""

        if not isinstance(label, AsymptoticSlowLabel):
            raise TypeError("label must be an AsymptoticSlowLabel")
        matches = tuple(
            binding
            for binding in self.physical_base_coverage.physical_base_bindings
            if binding.box_key == label.box_key
        )
        if len(matches) != 1:
            raise ValueError("signed label must resolve to exactly one physical-base binding")
        return matches[0]

    @property
    def covered_labels(self) -> frozenset[AsymptoticSlowLabel]:
        return self.family_estimate_coverage.covered_labels

    @property
    def covered_box_keys(self) -> frozenset[tuple[int, tuple[int, int, int]]]:
        return self.physical_base_coverage.covered_box_keys

    @property
    def source_key(self) -> tuple[str, str]:
        return self.family_estimate_coverage.scope.source_key

    def theorem_envelopes(self) -> dict[str, float]:
        """Return existing finite-scope theorem envelopes without re-fitting them."""

        return self.family_estimate_coverage.uniform_scope_bounds()

    @property
    def physical_estimate_scope_machine_checked(self) -> bool:
        """Mechanical coherence relative to the supplied theorem-facing scope."""

        return (
            self.physical_base_coverage.declared_scope_physical_base_coverage_machine_checked
            and self.family_estimate_coverage.scope_coverage_complete
            and self.physical_base_coverage.covered_signed_labels == self.covered_labels
            and self.physical_base_coverage.covered_box_keys
            == self.family_estimate_coverage.scope.box_keys
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
    def theorem_applications_machine_replayed(self) -> bool:
        return False

    @property
    def global_all_band_uniformity_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
