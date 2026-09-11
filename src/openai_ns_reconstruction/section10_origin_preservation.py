"""Section 10 late-origin localization identity, without a blow-up surrogate.

The pinned OpenAI formalization proves that the fixed spatial cutoff is one on
an open plateau containing the origin, that the time switch is one for
``t >= 3/4``, and consequently that the fully localized velocity agrees at the
origin with the incoming spatial curl on the late-time branch.  It then
transfers an *upstream proved* left-limit blow-up of that curl to the localized
velocity.

This module records the executable part of that theorem chain.  It certifies
only the exact localization factors at the origin for ``3/4 <= t < 1`` and the
associated cutoff product-rule algebra.  It does not sample a velocity, infer a
limit, construct the missing paper field, or certify that the upstream curl
actually diverges.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .coordinates import _finite
from .spatial_localization import (
    plateau_contains,
    section10_spatial_cutoff,
    section10_spatial_cutoff_gradient,
)
from .time_localization import LATE_START, time_switch
from .verify import finite_vector

BLOWUP_TIME = 1.0
ORIGIN = (0.0, 0.0, 0.0)


@dataclass(frozen=True)
class Section10LateOriginCertificate:
    """Fail-closed certificate for the Section 10 late-origin identity factors.

    The certificate is intentionally restricted to the pre-blow-up interval.
    The official pointwise identity itself is stated for all ``t >= 3/4``, but
    the transfer used for the singular limit approaches ``t=1`` from below.
    """

    t: float

    def __post_init__(self) -> None:
        t = _finite(self.t, "t")
        if t < LATE_START or not t < BLOWUP_TIME:
            raise ValueError("late-origin certificate requires 3/4 <= t < 1")
        object.__setattr__(self, "t", t)

        # These exact checks are regression guards.  If a future executable
        # cutoff no longer realizes the theorem's plateau geometry, the
        # certificate must fail rather than silently weakening the identity.
        if not plateau_contains(*ORIGIN):
            raise RuntimeError("origin is not in the Section 10 plateau")
        if self.spatial_cutoff_value != 1.0:
            raise RuntimeError("Section 10 spatial cutoff is not exactly one at the origin")
        if not np.array_equal(self.spatial_cutoff_gradient, np.zeros(3)):
            raise RuntimeError("Section 10 spatial cutoff gradient is not exactly zero at the origin")
        if self.time_switch_value != 1.0:
            raise RuntimeError("Section 10 time switch is not exactly one on the late branch")

    @property
    def origin_in_plateau(self) -> bool:
        return plateau_contains(*ORIGIN)

    @property
    def spatial_cutoff_value(self) -> float:
        return float(section10_spatial_cutoff(*ORIGIN, self.t))

    @property
    def spatial_cutoff_gradient(self) -> np.ndarray:
        return section10_spatial_cutoff_gradient(*ORIGIN, self.t)

    @property
    def time_switch_value(self) -> float:
        return float(time_switch(self.t))

    @property
    def identity_factors_exact(self) -> bool:
        """Whether the executable plateau/time factors are literally 1/0."""

        return (
            self.origin_in_plateau
            and self.spatial_cutoff_value == 1.0
            and np.array_equal(self.spatial_cutoff_gradient, np.zeros(3))
            and self.time_switch_value == 1.0
        )

    def cut_velocity_from_potential_data(
        self,
        spatial_curl_at_origin: object,
        potential_at_origin: object,
    ) -> np.ndarray:
        """Evaluate the cutoff product rule at the certified late origin.

        For supplied finite point data this evaluates

        ``c curl(A) + grad(c) x A``.

        The result must equal ``curl(A)`` because ``c=1`` and ``grad(c)=0`` at
        the origin.  Supplying point data is not a certificate that those data
        came from the paper's unresolved potential.
        """

        curl_value = finite_vector(spatial_curl_at_origin)
        potential_value = finite_vector(potential_at_origin)
        value = (
            self.spatial_cutoff_value * curl_value
            + np.cross(self.spatial_cutoff_gradient, potential_value)
        )
        return finite_vector(value)

    def activate_periodic_velocity(self, periodic_velocity_at_origin: object) -> np.ndarray:
        """Apply the official late time switch to a finite origin value.

        This isolates the time-localization part of the theorem chain.  The
        separate official no-overlap/plateau theorem is what identifies the
        periodic spatially localized value with the incoming curl.
        """

        value = finite_vector(periodic_velocity_at_origin)
        return finite_vector(self.time_switch_value * value)


def certify_section10_late_origin(t: float) -> Section10LateOriginCertificate:
    """Construct the fail-closed geometry/time certificate at one late time."""

    return Section10LateOriginCertificate(t)


def late_origin_cut_then_time_value(
    spatial_curl_at_origin: object,
    potential_at_origin: object,
    t: float,
) -> np.ndarray:
    """Compose executable cutoff-product and late-time identity factors.

    This is deliberately *not* named a blow-up verifier: a pointwise identity
    cannot establish the ``Tendsto ... atTop`` premise used by the paper/Lean
    blow-up theorem.  It is a reusable adapter for upstream certified data.
    """

    certificate = certify_section10_late_origin(t)
    cut_value = certificate.cut_velocity_from_potential_data(
        spatial_curl_at_origin, potential_at_origin
    )
    return finite_vector(certificate.time_switch_value * cut_value)
