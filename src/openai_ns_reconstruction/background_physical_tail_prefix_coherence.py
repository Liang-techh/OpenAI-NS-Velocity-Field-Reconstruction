"""Coherent increasing-truncation certificates for Section 5 physical tails.

The pinned all-order construction chooses one recursive SlowBorel/DiagonalScale
schedule and later selects request-dependent finite prefixes of that same
schedule.  A successful finite request therefore is not enough: when a stronger
physical jet/decay request forces a later truncation order, all already-retained
coefficient evidence and every already-selected cutoff scale must remain
literally unchanged.

This module checks that finite coherence property on top of the provider-owned
physical-tail certificates.  Both requests are constructed from the same real
hierarchy provider; callers do not supply coefficient rows, support tables,
recurrence identities, schedules, or omitted-order exponents.

For a weaker request (M0, P0) and a dominating request (M1, P1), M1 >= M0 and
P1 >= P0, it requires:

* one hierarchy/source/coefficient state;
* exact equality of all shared C[j,m] majorants, support witnesses, retained
  recurrence identities, and own-order dependency links;
* exact prefix equality of the recursive local scales and doubling-envelope
  scales; and
* the pinned first-omitted exponent gain h*(J1-J0) for every shared physical
  derivative row, together with the exact dyadic improvement 2^-(J1-J0).

This is still a two-request finite certificate.  It does not prove that the
provider is total, construct an infinite coefficient hierarchy, materialize the
PDE residual, or promote finite coherence to all-jets flatness/super-algebraic
convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .background_physical_tail_prefactor import (
    FinitePhysicalTailPrefactorCertificate,
    certify_finite_physical_tail_prefactor,
)
from .background_target_driven_exact_majorants import TargetDrivenExactMajorantProvider


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_SCALE_THEOREM = "NavierStokes.DiagonalScale.exists_diagonal_scales"
PINNED_ENVELOPE = "NavierStokes.DiagonalScale.doublingEnvelope"
PINNED_ADMISSIBLE_SCALES = "NavierStokes.SlowBorelBase.exists_admissibleScales"
PINNED_PHYSICAL_TAIL = "NavierStokes.SlowBorelBase.exists_physical_uncut_tail"


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


@dataclass(frozen=True)
class FiniteIncreasingPhysicalTailCertificate:
    """Two nested physical-tail requests sharing one literal finite schedule prefix."""

    earlier: FinitePhysicalTailPrefactorCertificate
    later: FinitePhysicalTailPrefactorCertificate
    formal_revision: str = PINNED_FORMAL_REVISION
    scale_theorem: str = PINNED_SCALE_THEOREM
    doubling_envelope: str = PINNED_ENVELOPE
    admissible_scales_theorem: str = PINNED_ADMISSIBLE_SCALES
    physical_tail_theorem: str = PINNED_PHYSICAL_TAIL

    def __post_init__(self) -> None:
        if not isinstance(self.earlier, FinitePhysicalTailPrefactorCertificate):
            raise TypeError("earlier must be a FinitePhysicalTailPrefactorCertificate")
        if not isinstance(self.later, FinitePhysicalTailPrefactorCertificate):
            raise TypeError("later must be a FinitePhysicalTailPrefactorCertificate")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.scale_theorem != PINNED_SCALE_THEOREM:
            raise ValueError("scale_theorem does not match exists_diagonal_scales")
        if self.doubling_envelope != PINNED_ENVELOPE:
            raise ValueError("doubling_envelope does not match the pinned recursive envelope")
        if self.admissible_scales_theorem != PINNED_ADMISSIBLE_SCALES:
            raise ValueError("admissible_scales_theorem does not match the pinned source")
        if self.physical_tail_theorem != PINNED_PHYSICAL_TAIL:
            raise ValueError("physical_tail_theorem does not match the pinned source")

        earlier_physical = self.earlier.first_omitted.uncut.physical
        later_physical = self.later.first_omitted.uncut.physical
        if self.earlier.first_omitted.h_exact != self.later.first_omitted.h_exact:
            raise ValueError("nested physical-tail requests use different h")
        if later_physical.max_derivative_order < earlier_physical.max_derivative_order:
            raise ValueError("later request must not reduce the physical derivative budget")
        if later_physical.target_physical_power < earlier_physical.target_physical_power:
            raise ValueError("later request must not reduce the target physical power")
        if later_physical.minimum_order != earlier_physical.minimum_order:
            raise ValueError("nested requests must use one common minimum truncation order")

        earlier_chain = self._exact_chain(self.earlier)
        later_chain = self._exact_chain(self.later)
        if earlier_chain.chain.hierarchy_id != later_chain.chain.hierarchy_id:
            raise ValueError("nested requests use different hierarchy ids")
        if earlier_chain.chain.source_revision != later_chain.chain.source_revision:
            raise ValueError("nested requests use different source revisions")
        if earlier_chain.chain.coefficient_state_id != later_chain.chain.coefficient_state_id:
            raise ValueError("nested requests use different coefficient states")

        J0 = self.earlier.truncation_order
        J1 = self.later.truncation_order
        if J1 < J0:
            raise ValueError("dominating physical request selected an earlier truncation order")

        earlier_prefix = earlier_chain.chain.supported_prefix
        later_prefix = later_chain.chain.supported_prefix
        if earlier_prefix.exact_support_box != later_prefix.exact_support_box:
            raise ValueError("nested requests use different common support boxes")

        # The exact-majorant tuple is ordered by (j,m), so the shorter request
        # must be a literal prefix of the stronger request.
        earlier_majorants = earlier_chain.exact_majorants
        if later_chain.exact_majorants[: len(earlier_majorants)] != earlier_majorants:
            raise ValueError("shared exact C[j,m] majorants changed across target extension")

        earlier_supports = earlier_prefix.supports.entries
        if later_prefix.supports.entries[: len(earlier_supports)] != earlier_supports:
            raise ValueError("shared coefficient support witnesses changed across target extension")

        earlier_links = earlier_chain.chain.links
        if later_chain.chain.links[: len(earlier_links)] != earlier_links:
            raise ValueError("shared own-order coefficient dependencies changed across target extension")

        earlier_cancellations = earlier_prefix.prefix.cancellations.identities
        later_cancellations = later_prefix.prefix.cancellations.identities
        if later_cancellations[: len(earlier_cancellations)] != earlier_cancellations:
            raise ValueError("shared retained recurrence identities changed across target extension")

        schedule0 = self.earlier.first_omitted.uncut.schedule
        schedule1 = self.later.first_omitted.uncut.schedule
        if schedule0.h != schedule1.h:
            raise ValueError("nested requests use different executable h")
        if schedule0.initial_lower_bound != schedule1.initial_lower_bound:
            raise ValueError("nested requests use different initial scale lower bounds")
        if schedule1.jet_bounds[: len(schedule0.jet_bounds)] != schedule0.jet_bounds:
            raise ValueError("recursive schedule C[j,m] rows are not prefix-coherent")
        if schedule1.local_scales[: len(schedule0.local_scales)] != schedule0.local_scales:
            raise ValueError("recursive local scales changed on a retained prefix")
        if schedule1.scales[: len(schedule0.scales)] != schedule0.scales:
            raise ValueError("doubling-envelope scales changed on a retained prefix")

        gain = self.physical_power_gain_exact
        earlier_rows = self.earlier.first_omitted.jet_powers
        later_rows = self.later.first_omitted.jet_powers
        shared = len(earlier_rows)
        if len(later_rows) < shared:
            raise RuntimeError("later request lost a previously requested derivative row")
        for row0, row1 in zip(earlier_rows, later_rows[:shared], strict=True):
            if row0.derivative_order != row1.derivative_order:
                raise RuntimeError("shared physical derivative rows are misaligned")
            if row1.physical_tail_power - row0.physical_tail_power != gain:
                raise ValueError("shared physical-tail exponent did not gain h*(J1-J0)")
            if row1.ordinary_tail_power - row0.ordinary_tail_power != gain:
                raise ValueError("shared ordinary-tail exponent did not gain h*(J1-J0)")

        expected_ratio = Fraction(1, 2) ** self.truncation_increment
        actual_ratio = (
            self.later.first_omitted.dyadic_prefactor_exact
            / self.earlier.first_omitted.dyadic_prefactor_exact
        )
        if actual_ratio != expected_ratio:
            raise ValueError("dyadic tail prefactor is not coherent under truncation extension")

        # A later retained stage cannot enlarge the exact uncut plateau.
        if self.later.first_omitted.uncut_radius_exact > self.earlier.first_omitted.uncut_radius_exact:
            raise ValueError("later recursive prefix unexpectedly enlarged the uncut radius")

    @staticmethod
    def _exact_chain(certificate: FinitePhysicalTailPrefactorCertificate):
        return certificate.first_omitted.uncut.physical.ordinary.chain

    @property
    def truncation_increment(self) -> int:
        return self.later.truncation_order - self.earlier.truncation_order

    @property
    def h_exact(self) -> Fraction:
        return self.earlier.first_omitted.h_exact

    @property
    def physical_power_gain_exact(self) -> Fraction:
        """Gain in every shared ordinary/physical exponent: ``h*(J1-J0)``."""
        return self.h_exact * self.truncation_increment

    @property
    def raw_first_omitted_power_gain_exact(self) -> Fraction:
        """Gain before diagonal half-power loss: ``2*h*(J1-J0)``."""
        return 2 * self.physical_power_gain_exact

    @property
    def dyadic_prefactor_ratio_exact(self) -> Fraction:
        """Return ``(2^-J1)/(2^-J0) = 2^-(J1-J0)`` exactly."""
        return Fraction(1, 2) ** self.truncation_increment

    @property
    def shared_schedule_prefix_verified(self) -> bool:
        return True

    @property
    def shared_coefficient_evidence_prefix_verified(self) -> bool:
        return True

    @property
    def increasing_truncation_exponent_gain_verified(self) -> bool:
        return True

    @property
    def finite_request_pair_only(self) -> bool:
        return True

    @property
    def provider_totality_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def actual_pde_residual_verified(self) -> bool:
        return False

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def all_jets_flat(self) -> bool:
        return False

    @property
    def super_algebraic(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False


def certify_increasing_physical_tail_pair(
    h: float,
    *,
    earlier_max_derivative_order: int,
    earlier_target_physical_power: Fraction | int,
    later_max_derivative_order: int,
    later_target_physical_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FiniteIncreasingPhysicalTailCertificate:
    """Build two nested provider-owned requests and verify literal prefix coherence."""

    earlier_max_derivative_order = _nonnegative_int(
        earlier_max_derivative_order, "earlier_max_derivative_order"
    )
    later_max_derivative_order = _nonnegative_int(
        later_max_derivative_order, "later_max_derivative_order"
    )
    if later_max_derivative_order < earlier_max_derivative_order:
        raise ValueError("later_max_derivative_order must be >= earlier_max_derivative_order")

    earlier = certify_finite_physical_tail_prefactor(
        h,
        max_derivative_order=earlier_max_derivative_order,
        target_physical_power=earlier_target_physical_power,
        provider=provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    later = certify_finite_physical_tail_prefactor(
        h,
        max_derivative_order=later_max_derivative_order,
        target_physical_power=later_target_physical_power,
        provider=provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    return FiniteIncreasingPhysicalTailCertificate(earlier=earlier, later=later)
