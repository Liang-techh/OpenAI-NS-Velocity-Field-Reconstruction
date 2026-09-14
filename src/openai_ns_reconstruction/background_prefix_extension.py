"""Successive finite-prefix closure for Section 5 SlowBorel residual bounds.

The repository already has three separate fail-closed ingredients:

* dependency-backed finite-prefix ``C[j,m]`` bounds feeding the pinned
  SlowBorel/DiagonalScale recursive cutoff schedule;
* exact retained-recurrence cancellation certificates, with no numerical
  near-zero gate; and
* the exact first-omitted-slow-order residual majorant.

This module couples those ingredients for one *successive* truncation step
``N -> N+1`` and for an arbitrary chosen *finite* chain of such steps.  A step
is accepted only when the hierarchy provenance, old cutoff rows/scales, and old
exact recurrence identities are literally preserved, the first omitted exponent
gains exactly ``2*h``, and the certified finite residual majorant at the same
``q`` does not increase.  A finite chain additionally requires exact equality
of every shared cutoff/residual endpoint, so independently valid but cross-wired
prefixes cannot be spliced together.

The resulting objects are still finite-prefix ``formal-structure``.  They do not
manufacture hierarchy bounds or recurrence identities.  Quantifying over an
arbitrary user-supplied finite chain is not an infinite coherent family and is
not by itself Proposition 5.3, all-jets-flatness, or super-algebraic convergence.
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


@dataclass(frozen=True)
class FiniteCoherentSlowBorelResidualChain:
    """Compose any nonempty *finite* chain of certified one-order extensions.

    Each step is already fail-closed.  This class adds the global coherence
    requirement that the extended endpoint of step ``k`` is literally the
    previous endpoint of step ``k+1`` for both the hierarchy-backed cutoff
    witness and the exact-cancellation residual witness.  It also rechecks a
    single hierarchy/source/coefficient-state identity, one physical
    ``q/h/base_power`` regime, contiguous truncation orders, the exact cumulative
    first-omitted exponent gain, and endpoint majorant monotonicity.

    The quantifier here is only over the finite tuple supplied by the caller.
    No property of this object certifies extendability for every order.
    """

    steps: tuple[SuccessiveSlowBorelResidualExtension, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.steps, tuple):
            raise TypeError("steps must be a tuple of certified successive extensions")
        if not self.steps:
            raise ValueError("finite coherent SlowBorel chain must contain at least one step")
        if any(not isinstance(step, SuccessiveSlowBorelResidualExtension) for step in self.steps):
            raise TypeError("every chain step must be a SuccessiveSlowBorelResidualExtension")

        first = self.steps[0]
        last = self.steps[-1]
        first_bounds = first.previous_cutoff.bounds
        first_cancellations = first.previous_residual.cancellations
        first_tail = first.previous_residual.tail

        for index, step in enumerate(self.steps):
            bounds = step.previous_cutoff.bounds
            cancellations = step.previous_residual.cancellations
            tail = step.previous_residual.tail
            if bounds.hierarchy_id != first_bounds.hierarchy_id:
                raise ValueError("finite chain crosses hierarchy identities")
            if bounds.source_revision != first_bounds.source_revision:
                raise ValueError("finite chain crosses source revisions")
            if cancellations.coefficient_state_id != first_cancellations.coefficient_state_id:
                raise ValueError("finite chain crosses coefficient states")
            if tail.q != first_tail.q:
                raise ValueError("finite chain changes q")
            if tail.h != first_tail.h:
                raise ValueError("finite chain changes h")
            if tail.base_power != first_tail.base_power:
                raise ValueError("finite chain changes base power")

            if index:
                previous = self.steps[index - 1]
                if previous.extended_cutoff != step.previous_cutoff:
                    raise ValueError("adjacent steps do not share the exact cutoff endpoint")
                if previous.extended_residual != step.previous_residual:
                    raise ValueError("adjacent steps do not share the exact residual endpoint")
                if previous.extended_order != step.previous_order:
                    raise ValueError("adjacent steps are not contiguous in truncation order")

        if self.end_order - self.start_order != self.step_count:
            raise ValueError("finite chain does not advance exactly one order per step")

        gain = (
            last.extended_residual.tail.exponent_exact
            - first.previous_residual.tail.exponent_exact
        )
        expected_gain = slow_order_exact(first_tail.h, self.step_count)
        if gain != expected_gain:
            raise ValueError("finite chain first-omitted exponent gain is not exactly 2*h per step")

        if math.isnan(self.initial_log_majorant) or math.isnan(self.final_log_majorant):
            raise ValueError("finite-chain endpoint residual log-majorants must not be NaN")
        if self.final_log_majorant > self.initial_log_majorant:
            raise ValueError("finite-chain endpoint residual majorant increased")

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def infinite_coherent_family(self) -> bool:
        return False

    @property
    def step_count(self) -> int:
        return len(self.steps)

    @property
    def start_order(self) -> int:
        return self.steps[0].previous_order

    @property
    def end_order(self) -> int:
        return self.steps[-1].extended_order

    @property
    def initial_log_majorant(self) -> float:
        return self.steps[0].previous_residual.tail.log_bound

    @property
    def final_log_majorant(self) -> float:
        return self.steps[-1].extended_residual.tail.log_bound

    @property
    def total_first_omitted_exponent_gain(self) -> Fraction:
        return (
            self.steps[-1].extended_residual.tail.exponent_exact
            - self.steps[0].previous_residual.tail.exponent_exact
        )

    @property
    def strict_improvement_count(self) -> int:
        return sum(step.majorant_strictly_improves for step in self.steps)

    @property
    def overall_majorant_improvement(self) -> float:
        return self.initial_log_majorant - self.final_log_majorant


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


def certify_finite_coherent_slow_borel_residual_chain(
    steps: tuple[SuccessiveSlowBorelResidualExtension, ...],
) -> FiniteCoherentSlowBorelResidualChain:
    """Compose certified successive steps without making an infinite-order claim."""

    return FiniteCoherentSlowBorelResidualChain(steps=steps)
