"""Physical late-time slab bridge for Section 9 local-sum admission.

This module connects the exact-rational Eq. (9.21) q-strip admission on main to
one *finite* physical-time slab inside the fixed Section 10 support cylinder.
The bridge uses only Eq. (4.1) and the official support bound ``|z|<=1/4``.

For ``tau=1-t`` and ``0<h<1/2`` the physical root satisfies

    tau = q - z^2 q^(2h).

On the Section 10 late plateau ``t>=3/4`` and support cylinder ``z^2<=1/16``
we first have ``q<1``: otherwise ``q^(2h)<=q`` would imply
``tau >= q(1-z^2) >= 15/16``, contradicting ``tau<=1/4``.  Hence
``q^(2h)<=1`` and therefore

    1-t <= q <= (1-t) + z^2 <= (1-t) + 1/16.

Thus on a closed pre-endpoint slab ``t0<=t<=t1<1`` the whole official spatial
support lies in the exact q-strip

    1-t1 <= q <= 1-t0+1/16.

If that upper bound is strictly below the independently certified common
``q_big`` from :mod:`section9_local_sum_admission`, the landed local-finiteness
arithmetic applies uniformly on this physical slab.  The construction is
intentionally fail-closed at ``t1=1``: the lower q bound then collapses to zero,
so this finite-slab bridge is *not* an endpoint convergence theorem and cannot
manufacture the Section 10 all-order residual majorants.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Rational, Real
import math

from .section9_local_sum_admission import (
    Section9LocalSumAdmissionCertificate,
    Section9LocalSumHypothesis,
    certify_eq_9_21_local_finiteness,
)
from .spatial_localization import SUPPORT_HALF_HEIGHT


LATE_START_EXACT = Fraction(3, 4)
SECTION10_ENDPOINT_EXACT = Fraction(1, 1)
SUPPORT_HALF_HEIGHT_EXACT = Fraction(1, 4)
SUPPORT_Z_SQUARED_EXACT = Fraction(1, 16)


def _finite_fraction(value: object, name: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real value")
    if isinstance(value, Fraction):
        out = value
    elif isinstance(value, Integral):
        out = Fraction(int(value), 1)
    elif isinstance(value, Rational):
        out = Fraction(value)
    elif isinstance(value, Real):
        x = float(value)
        if not math.isfinite(x):
            raise ValueError(f"{name} must be finite")
        out = Fraction.from_float(x)
    else:
        raise TypeError(f"{name} must be a real rationalizable value")
    return out


@dataclass(frozen=True)
class Section9PhysicalSlabLocalSumCertificate:
    """Uniform Eq. (9.21) local-finiteness admission on one physical slab.

    The nested ``local_sum`` certificate is the exact-rational Stage-9
    arithmetic already present on main.  This wrapper additionally certifies
    that the *entire fixed Section 10 support cylinder* maps into that q-strip
    for ``t_lower<=t<=t_upper<1`` by the analytic Eq. (4.1) bounds above.

    No actual correction field is constructed here, and the endpoint is
    deliberately excluded.
    """

    t_lower: Fraction
    t_upper: Fraction
    q_lower_bound: Fraction
    q_upper_bound: Fraction
    local_sum: Section9LocalSumAdmissionCertificate
    support_geometry_verified: bool = True
    similarity_bound_verified: bool = True
    endpoint_covered: bool = False
    actual_correction_fields_verified: bool = False
    eq_9_21_sum_constructed: bool = False
    endpoint_uniform_residual_majorants_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def max_potentially_active_stage(self) -> int:
        return self.local_sum.max_potentially_active_stage

    @property
    def first_guaranteed_inactive_stage(self) -> int:
        return self.local_sum.first_guaranteed_inactive_stage


def certify_section10_physical_slab_local_finiteness(
    hypothesis: Section9LocalSumHypothesis,
    *,
    t_lower: Fraction | float | int,
    t_upper: Fraction | float | int,
) -> Section9PhysicalSlabLocalSumCertificate:
    """Transfer the Section 9 q-strip admission to the fixed Section 10 support.

    ``t_lower`` and ``t_upper`` describe a closed physical-time slab.  The
    lower endpoint must already lie in the official ``t>=3/4`` unit plateau,
    while the upper endpoint must stay strictly below the singular time 1.

    The bridge refuses a slab unless the analytic q upper bound for *every*
    point in the official support cylinder lies strictly inside the supplied
    common q-domain.  It never estimates that domain from samples.
    """
    if not isinstance(hypothesis, Section9LocalSumHypothesis):
        raise TypeError("hypothesis must be a Section9LocalSumHypothesis")

    lower = _finite_fraction(t_lower, "t_lower")
    upper = _finite_fraction(t_upper, "t_upper")
    if lower < LATE_START_EXACT:
        raise ValueError("physical slab must start in the official t>=3/4 plateau")
    if not lower < upper:
        raise ValueError("physical slab requires t_lower < t_upper")
    if upper >= SECTION10_ENDPOINT_EXACT:
        raise ValueError("physical slab must end strictly before t=1; endpoint control is separate")

    if SUPPORT_HALF_HEIGHT != float(SUPPORT_HALF_HEIGHT_EXACT):
        raise ArithmeticError("Section 10 support geometry drifted from |z|<=1/4")

    q_lower = SECTION10_ENDPOINT_EXACT - upper
    q_upper = SECTION10_ENDPOINT_EXACT - lower + SUPPORT_Z_SQUARED_EXACT
    if not q_lower > 0:
        raise ArithmeticError("pre-endpoint slab must give a positive uniform q lower bound")
    if not q_upper < hypothesis.q_big:
        raise ValueError(
            "analytic q upper bound for the fixed Section 10 support is not inside q_big"
        )

    local_sum = certify_eq_9_21_local_finiteness(
        hypothesis,
        q_lower=q_lower,
        q_upper=q_upper,
    )
    return Section9PhysicalSlabLocalSumCertificate(
        t_lower=lower,
        t_upper=upper,
        q_lower_bound=q_lower,
        q_upper_bound=q_upper,
        local_sum=local_sum,
    )
