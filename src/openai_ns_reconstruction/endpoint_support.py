"""Pointwise compact-support certificates for the Section 10 endpoint glue.

The pinned ``SpatialBorelExtension.extension_zero_of_coefficients_zero`` theorem
says that a spatial point at which every Taylor--Borel coefficient vanishes
stays identically zero under the extension.  Executably, at any fixed future
time ``t>T`` only finitely many coefficients can contribute because the landed
``DoublingEnvelope`` eventually moves every later cutoff outside its support.

This module records that finite pointwise consequence for the repository's
formal endpoint-force path.  It does not prove that the actual residual jets
vanish outside the Section 10 cylinder; that remains an upstream smoothness and
support obligation.  It also does not turn caller-supplied jets into a
paper-exact force.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np

from .coordinates import _finite
from .cutoffs import standard_cutoff
from .endpoint_borel import BorelRightExtension, LocalScale, borel_monomial_coefficient
from .endpoint_jets import FullSpacetimeJetFamily, normal_jet_family
from .spatial_localization import support_cylinder_contains
from .verify import finite_vector

VectorField = Callable[[float, float, float, float], np.ndarray]


def _zero_vector(value: object) -> bool:
    vector = finite_vector(value)
    return bool(np.array_equal(vector, np.zeros(3)))


def _scaled_argument(scale: int, s: float) -> float:
    """Mirror the Borel evaluator's support comparison without clipping scales."""
    try:
        return float(scale) * s
    except OverflowError:
        return math.inf


@dataclass(frozen=True)
class FormalForceSupportPointCertificate:
    """Finite certificate that one formal-force value is zero outside support.

    ``checked_degrees`` lists precisely the normal Taylor coefficients whose
    cutoff and monomial factors can be nonzero at the queried future point.  On
    the past branch it is empty because the supplied past residual is checked
    directly.  ``certified`` is deliberately pointwise; it is not a proof of a
    global compact-support statement over all spacetime points.
    """

    x: float
    y: float
    z: float
    t: float
    endpoint: float
    branch: str
    checked_degrees: tuple[int, ...]
    relevant_inputs_zero: bool
    formal_value_zero: bool

    @property
    def certified(self) -> bool:
        branch_time_ok = (
            (self.branch == "past" and self.t < self.endpoint)
            or (self.branch == "endpoint" and self.t == self.endpoint)
            or (self.branch == "future" and self.t > self.endpoint)
        )
        return (
            branch_time_ok
            and not support_cylinder_contains(self.x, self.y, self.z)
            and self.relevant_inputs_zero
            and self.formal_value_zero
        )


def _future_certificate(
    full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    endpoint: float,
) -> FormalForceSupportPointCertificate:
    normal = normal_jet_family(full_jets)
    right = BorelRightExtension(normal, local_scale, endpoint=endpoint)

    if t == endpoint:
        coefficient = normal(0, x, y, z)
        value = right(x, y, z, t)
        return FormalForceSupportPointCertificate(
            x=x,
            y=y,
            z=z,
            t=t,
            endpoint=endpoint,
            branch="endpoint",
            checked_degrees=(0,),
            relevant_inputs_zero=_zero_vector(coefficient),
            formal_value_zero=_zero_vector(value),
        )

    s = t - endpoint
    checked: list[int] = []
    all_zero = True
    degree = 0
    while degree < 4096:
        scale = right.scale(degree)
        argument = _scaled_argument(scale, s)
        if argument >= 1.0:
            break
        cutoff = standard_cutoff(argument)
        monomial = borel_monomial_coefficient(s, degree)
        if cutoff != 0.0 and monomial != 0.0:
            coefficient = normal(degree, x, y, z)
            checked.append(degree)
            all_zero = all_zero and _zero_vector(coefficient)
        degree += 1
    else:
        raise RuntimeError("doubling-envelope local-finiteness guard was exceeded")

    value = right(x, y, z, t)
    return FormalForceSupportPointCertificate(
        x=x,
        y=y,
        z=z,
        t=t,
        endpoint=endpoint,
        branch="future",
        checked_degrees=tuple(checked),
        relevant_inputs_zero=all_zero,
        formal_value_zero=_zero_vector(value),
    )


def formal_force_support_point_certificate(
    past_residual: VectorField,
    full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    endpoint: float = 1.0,
) -> FormalForceSupportPointCertificate:
    """Certify one outside-cylinder value of the formal Section 10 force.

    The branch convention matches ``traced_residual.formal_force_from_full_spacetime_jets``:

    * ``t<T`` checks the supplied past residual directly;
    * ``t=T`` checks only the degree-zero endpoint coefficient; and
    * ``t>T`` checks exactly the finitely many normal coefficients that can
      contribute to the Taylor--Borel sum at that time.

    The queried spatial point must lie outside the official closed Section 10
    support cylinder.  Exact zeros are required: no tolerance, sampling-based
    extrapolation, or fitted support boundary is accepted.
    """
    if not callable(past_residual):
        raise ValueError("past_residual must be callable")
    if not callable(full_jets):
        raise ValueError("full_jets must be callable")
    if not callable(local_scale):
        raise ValueError("local_scale must be callable")

    x = _finite(x, "x")
    y = _finite(y, "y")
    z = _finite(z, "z")
    t = _finite(t, "t")
    endpoint = _finite(endpoint, "endpoint")
    if support_cylinder_contains(x, y, z):
        raise ValueError("support certificate requires a point outside the closed Section 10 cylinder")

    if t < endpoint:
        value = finite_vector(past_residual(x, y, z, t))
        is_zero = _zero_vector(value)
        return FormalForceSupportPointCertificate(
            x=x,
            y=y,
            z=z,
            t=t,
            endpoint=endpoint,
            branch="past",
            checked_degrees=(),
            relevant_inputs_zero=is_zero,
            formal_value_zero=is_zero,
        )

    return _future_certificate(
        full_jets,
        local_scale,
        x,
        y,
        z,
        t,
        endpoint=endpoint,
    )
