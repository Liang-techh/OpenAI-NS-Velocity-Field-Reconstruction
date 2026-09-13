"""Successive finite-prefix closure for Section 5 SlowBorel residual bounds.

The repository already has three separate fail-closed ingredients:

* dependency-backed finite-prefix ``C[j,m]`` bounds feeding the pinned
  SlowBorel/DiagonalScale recursive cutoff schedule;
* exact retained-recurrence cancellation certificates, with no numerical
  near-zero gate; and
* the exact first-omitted-slow-order residual majorant.

This module couples those ingredients for one *successive* truncation step
``N -> N+1``.  It accepts an extension only when the hierarchy provenance,
old cutoff rows/scales, and old exact recurrence identities are literally
preserved, the first omitted exponent gains exactly ``2*h``, and the certified
finite residual majorant at the same ``q`` does not increase.

The resulting object is still finite-prefix ``formal-structure``.  It does not
manufacture hierarchy bounds or recurrence identities, and a sequence of such
finite witnesses is not by itself Proposition 5.3, all-jets-flatness, or
super-algebraic convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from .background_hierarchy_cutoff_bounds import HierarchySlowBorelPrefixWitness
from .background_recurrence_cancellation import CertifiedFullResidualTailMajorant
from .background_truncation_residual import slow_order_exact


@dataclass(frozen=True)
class SuccessiveSlowBorelResidualExtension:
    """Machine-checkable one-order extension ``N -> N+1``.

    All four inputs must already have passed their own fail-closed constructors.
    This class checks the cross-layer facts that matter for increasing
    truncation:

    * the extended ``C[j,m]`` certificate is an exact prefix extension of the
      previous one and therefore owns the same old SlowBorel scales;
    * the extended exact-cancellation certificate preserves every old retained
      recurrence identity verbatim;
    * cutoff and recurrence evidence refer to one hierarchy/source revision;
    * both tail majorants use the same ``q``, ``h`` and physical base power;
    * the first omitted slow exponent increases by exactly ``2*h``; and
    * the finite certified full-residual log-majorant does not increase.

    The last condition is intentionally checked rather than inferred from the
    exponent alone, because the omitted coefficient prefactor may grow with
    truncation order.
    """

    previous_cutoff: HierarchySlowBorelPrefixWitness
    extended_cutoff: HierarchySlowBorelPrefixWitness
    previous_residual: CertifiedFullResidualTailMajorant
    extended_residual: CertifiedFullResidualTailMajorant

    def __post_init__(self) -> None:
        if not isinstance(self.previous_cutoff, HierarchySlowBorelPrefixWitness):
            raise TypeError("previous_cutoff must be a HierarchySlowBorelPrefixWitness")
        if not isinstance(self.extended_cutoff, HierarchySlowBorelPrefixWitness):
            raise TypeError("extended_cutoff must be a HierarchySlowBorelPrefixWitness")
        if not isinstance(self.previous_residual, CertifiedFullResidualTailMajorant):
            raise TypeError("previous_residual must be a CertifiedFullResidualTailMajorant")
        if not isinstance(self.extended_residual, CertifiedFullResidualTailMajorant):
            raise TypeError("extended_residual must be a CertifiedFullResidualTailMajorant")

        previous_order = self.previous_residual.cancellations.max_order
        extended_order = self.extended_residual.cancellations.max_order
        if previous_order < 1:
            raise ValueError("successive SlowBorel extension currently requires N >= 1")
        if extended_order != previous_order + 1:
            raise ValueError("residual certificate must extend the truncation by exactly one order")
        if self.previous_cutoff.bounds.max_order != previous_order:
            raise ValueError("previous cutoff order does not match previous cancellation order")
        if self.extended_cutoff.bounds.max_order != extended_order:
            raise ValueError("extended cutoff order does not match extended cancellation order")

        old_bounds = self.previous_cutoff.bounds
        new_bounds = self.extended_cutoff.bounds
        old_cancellations = self.previous_residual.cancellations
        new_cancellations = self.extended_residual.cancellations

        if old_bounds.hierarchy_id != new_bounds.hierarchy_id:
            raise ValueError("cutoff extension crosses hierarchy identities")
        if old_bounds.source_revision != new_bounds.source_revision:
            raise ValueError("cutoff extension crosses source revisions")
        if old_cancellations.hierarchy_id != new_cancellations.hierarchy_id:
            raise ValueError("recurrence extension crosses hierarchy identities")
        if old_cancellations.source_revision != new_cancellations.source_revision:
            raise ValueError("recurrence extension crosses source revisions")
        if old_cancellations.coefficient_state_id != new_cancellations.coefficient_state_id:
            raise ValueError("recurrence extension crosses coefficient states")
        if old_bounds.hierarchy_id != old_cancellations.hierarchy_id:
            raise ValueError("cutoff bounds and recurrence identities use different hierarchies")
        if old_bounds.source_revision != old_cancellations.source_revision:
            raise ValueError("cutoff bounds and recurrence identities use different source revisions")

        prefix_bound_count = len(old_bounds.entries)
        if new_bounds.entries[:prefix_bound_count] != old_bounds.entries:
            raise ValueError("extended C[j,m] certificate changes a previously certified bound")
        if new_cancellations.identities[: previous_order + 1] != old_cancellations.identities:
            raise ValueError("extended recurrence certificate changes a retained identity")

        old_schedule = self.previous_cutoff.schedule
        new_schedule = self.extended_cutoff.schedule
        if old_schedule.h != new_schedule.h:
            raise ValueError("cutoff extension changes h")
        if old_schedule.initial_lower_bound != new_schedule.initial_lower_bound:
            raise ValueError("cutoff extension changes the stage-zero lower bound")
        if new_schedule.jet_bounds[:previous_order] != old_schedule.jet_bounds:
            raise ValueError("extended scheduler changes a previously admitted C[j,m] row")
        if new_schedule.local_scales[: previous_order + 1] != old_schedule.local_scales:
            raise ValueError("extended scheduler changes a previous local scale")
        if new_schedule.scales[: previous_order + 1] != old_schedule.scales:
            raise ValueError("extended scheduler changes a previous recursive scale")

        old_tail = self.previous_residual.tail
        new_tail = self.extended_residual.tail
        if old_tail.order != previous_order or new_tail.order != extended_order:
            raise ValueError("tail majorant order does not match its cancellation certificate")
        if old_tail.q != new_tail.q:
            raise ValueError("successive residual comparison requires the same q")
        if old_tail.h != new_tail.h or old_tail.h != old_schedule.h:
            raise ValueError("cutoff and residual witnesses must use one common h")
        if old_tail.base_power != new_tail.base_power:
            raise ValueError("successive residual comparison requires the same base power")

        gain = new_tail.exponent_exact - old_tail.exponent_exact
        expected_gain = slow_order_exact(old_tail.h, 1)
        if gain != expected_gain:
            raise ValueError("first omitted exponent did not gain exactly 2*h")

        if math.isnan(old_tail.log_bound) or math.isnan(new_tail.log_bound):
            raise ValueError("residual log-majorants must not be NaN")
        if new_tail.log_bound > old_tail.log_bound:
            raise ValueError(
                "extended finite residual majorant is larger at the certified q; "
                "exponent gain alone cannot certify truncation improvement"
            )

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def previous_order(self) -> int:
        return self.previous_residual.cancellations.max_order

    @property
    def extended_order(self) -> int:
        return self.extended_residual.cancellations.max_order

    @property
    def previous_first_omitted_order(self) -> int:
        return self.previous_residual.first_omitted_order

    @property
    def extended_first_omitted_order(self) -> int:
        return self.extended_residual.first_omitted_order

    @property
    def exponent_gain_exact(self) -> Fraction:
        return (
            self.extended_residual.tail.exponent_exact
            - self.previous_residual.tail.exponent_exact
        )

    @property
    def new_recursive_scale(self) -> int:
        return self.extended_cutoff.schedule.scales[-1]

    @property
    def majorant_nonincreasing(self) -> bool:
        return self.extended_residual.tail.log_bound <= self.previous_residual.tail.log_bound

    @property
    def majorant_strictly_improves(self) -> bool:
        return self.extended_residual.tail.log_bound < self.previous_residual.tail.log_bound


def certify_successive_slow_borel_residual_extension(
    previous_cutoff: HierarchySlowBorelPrefixWitness,
    extended_cutoff: HierarchySlowBorelPrefixWitness,
    previous_residual: CertifiedFullResidualTailMajorant,
    extended_residual: CertifiedFullResidualTailMajorant,
) -> SuccessiveSlowBorelResidualExtension:
    """Certify one provenance-preserving, exact-cancellation truncation step."""

    return SuccessiveSlowBorelResidualExtension(
        previous_cutoff=previous_cutoff,
        extended_cutoff=extended_cutoff,
        previous_residual=previous_residual,
        extended_residual=extended_residual,
    )
