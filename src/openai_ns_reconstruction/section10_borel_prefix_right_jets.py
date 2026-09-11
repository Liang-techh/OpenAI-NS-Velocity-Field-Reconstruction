"""Finite-prefix endpoint-jet algebra for the Section 10 Borel future branch.

The pinned ``SpatialBorelExtension.rightExtension_right_jets`` theorem says that
all one-sided time derivatives of the *infinite* Taylor--Borel right extension
at ``t=T`` equal the prescribed normal jets.  The full theorem also uses the
joint-smoothness / summable-derivative machinery in ``SpatialBorelExtension``.

This module implements one strictly smaller, reusable piece of that argument.
For degrees ``0..N`` it records the exact common interval on which every cutoff
appearing in the finite prefix is on its unit plateau.  On that interval the
prefix is exactly

    sum_{j=0}^N (t-T)^j / j! * a_j(x),

so its degree-``n`` endpoint time derivative is exactly ``a_n(x)``.  Here
``a_n`` is not an arbitrary vector: it is the all-time-direction contraction of
the supplied dense full spacetime endpoint tensor, using the repository's
``(time,x,y,z)`` convention.

This finite-prefix statement is intentionally *not* promoted to the infinite
Borel theorem.  In particular, caller-supplied endpoint tensors and scales do
not prove that the tensors are locally-uniform limits of the actual Section 9
residual, that the official analytic template bounds hold, or that derivative
tails may be interchanged with the infinite sum.  ``paper_exact_velocity``
therefore remains unavailable.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import operator
from typing import Iterable

import numpy as np

from .coordinates import _finite
from .endpoint_borel import BorelRightExtension, LocalScale
from .endpoint_jets import FullSpacetimeJetFamily, normal_jet_family, normal_time_jet_from_full

SECTION10_ENDPOINT = 1.0


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(out)


def _readonly_vector(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    out = np.array(out, dtype=float, copy=True)
    out.setflags(write=False)
    return out


def _prefix_endpoint_derivative(
    normal_jets: tuple[np.ndarray, ...], degree: int
) -> np.ndarray:
    """Differentiate the degree-``degree`` Taylor prefix at offset zero.

    This is deliberately written as prefix algebra rather than simply returning
    ``normal_jets[degree]``.  For ``j<=degree``, the degree-th derivative at
    zero of ``s^j/j!`` is zero unless ``j=degree``, when it is one.
    """

    out = np.zeros(3, dtype=float)
    for j, jet in enumerate(normal_jets[: degree + 1]):
        if j == degree:
            out += jet
        # Every j<degree term has an exact zero degree-th derivative at s=0.
    return _readonly_vector(out, "finite-prefix endpoint time jet")


@dataclass(frozen=True)
class BorelPrefixRightJetDegreeCertificate:
    """Exact finite-prefix right-jet statement for one derivative degree."""

    degree: int
    point: tuple[float, float, float]
    endpoint: float
    prefix_scales: tuple[int, ...]
    common_plateau_radius: Fraction
    prescribed_normal_jet: np.ndarray
    prefix_endpoint_time_jet: np.ndarray

    def __post_init__(self) -> None:
        degree = _natural(self.degree, "degree")
        endpoint = _finite(self.endpoint, "endpoint")
        point = tuple(_finite(value, name) for value, name in zip(self.point, ("x", "y", "z")))
        scales = tuple(int(scale) for scale in self.prefix_scales)
        if len(scales) != degree + 1 or any(scale < 1 for scale in scales):
            raise ValueError("prefix_scales must contain one positive scale for every degree 0..n")
        if any(scales[j + 1] < 2 * scales[j] for j in range(len(scales) - 1)):
            raise ValueError("prefix scales must obey the landed doubling envelope")
        radius = self.common_plateau_radius
        if not isinstance(radius, Fraction) or radius <= 0:
            raise ValueError("common_plateau_radius must be a positive exact Fraction")
        expected = Fraction(1, 2 * max(scales))
        if radius != expected:
            raise ValueError("common plateau radius must be exactly 1/(2 max prefix scale)")

        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "point", point)
        object.__setattr__(self, "endpoint", endpoint)
        object.__setattr__(self, "prefix_scales", scales)
        object.__setattr__(self, "prescribed_normal_jet", _readonly_vector(self.prescribed_normal_jet, "prescribed normal jet"))
        object.__setattr__(self, "prefix_endpoint_time_jet", _readonly_vector(self.prefix_endpoint_time_jet, "prefix endpoint time jet"))

    @property
    def unit_plateau_verified(self) -> bool:
        """All prefix cutoffs are one for ``0 <= t-T <= radius``.

        The pinned/executable cutoff is exactly one for absolute argument at
        most ``1/2``.  The check is exact rational arithmetic and never rounds a
        tiny plateau radius to binary64 zero.
        """

        return all(
            Fraction(scale, 1) * self.common_plateau_radius <= Fraction(1, 2)
            for scale in self.prefix_scales
        )

    @property
    def finite_prefix_right_jet_verified(self) -> bool:
        return (
            self.endpoint == SECTION10_ENDPOINT
            and self.unit_plateau_verified
            and np.array_equal(self.prefix_endpoint_time_jet, self.prescribed_normal_jet)
        )


@dataclass(frozen=True)
class Section10BorelPrefixRightJetsCertificate:
    """Finite family of exact prefix endpoint-normal-jet identities.

    The negative flags are part of the contract: this object is allowed to
    certify finite Taylor-prefix algebra only.  It cannot be interpreted as the
    infinite ``rightExtension_right_jets`` theorem or as a smooth force.
    """

    point: tuple[float, float, float]
    endpoint: float
    max_degree: int
    rows: tuple[BorelPrefixRightJetDegreeCertificate, ...]
    dense_full_jet_shape_verified: bool = True
    finite_prefix_right_jet_algebra_verified: bool = True
    actual_section9_residual_limits_verified: bool = False
    analytic_template_bounds_verified: bool = False
    infinite_borel_right_jets_verified: bool = False
    all_order_borel_smoothness_verified: bool = False
    smooth_compact_forcing_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_prefix_ready(self) -> bool:
        expected_degrees = tuple(range(self.max_degree + 1))
        return (
            self.endpoint == SECTION10_ENDPOINT
            and tuple(row.degree for row in self.rows) == expected_degrees
            and all(row.point == self.point and row.endpoint == self.endpoint for row in self.rows)
            and all(row.finite_prefix_right_jet_verified for row in self.rows)
            and self.dense_full_jet_shape_verified
            and self.finite_prefix_right_jet_algebra_verified
            and not self.actual_section9_residual_limits_verified
            and not self.analytic_template_bounds_verified
            and not self.infinite_borel_right_jets_verified
            and not self.all_order_borel_smoothness_verified
            and not self.smooth_compact_forcing_verified
            and not self.paper_exact_velocity_available
        )


def certify_section10_borel_prefix_right_jets(
    full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    x: object,
    y: object,
    z: object,
    max_degree: object,
    *,
    endpoint: object = SECTION10_ENDPOINT,
) -> Section10BorelPrefixRightJetsCertificate:
    """Certify the exact finite-prefix right jets through ``max_degree``.

    ``local_scale`` is passed through the same landed ``DoublingEnvelope`` used
    by :class:`~openai_ns_reconstruction.endpoint_borel.BorelRightExtension`.
    For degree ``n`` the common plateau of terms ``0..n`` has exact radius
    ``1/(2*max_{j<=n} a_j)``.  Within that right neighborhood all those cutoff
    factors are identically one, so finite Taylor algebra gives the endpoint
    time jet exactly.

    The function validates every dense full spacetime tensor that it consumes.
    It does not sample the future branch, fit endpoint derivatives, or infer any
    missing analytic bounds from numerical values.
    """

    if not callable(full_jets):
        raise ValueError("full_jets must be callable")
    if not callable(local_scale):
        raise ValueError("local_scale must be callable")
    degree_max = _natural(max_degree, "max_degree")
    endpoint_value = _finite(endpoint, "endpoint")
    if endpoint_value != SECTION10_ENDPOINT:
        raise ValueError("Section 10 Borel right-jet certificate requires endpoint=1")
    point = (_finite(x, "x"), _finite(y, "y"), _finite(z, "z"))

    # Validate and contract each full tensor exactly once.  These are the
    # prescribed normal coefficients used by the pinned SpacetimeGluing path.
    normal_jets = tuple(
        _readonly_vector(
            normal_time_jet_from_full(full_jets, degree, *point),
            f"normal jet {degree}",
        )
        for degree in range(degree_max + 1)
    )

    right = BorelRightExtension(
        normal_jet_family(full_jets), local_scale, endpoint=endpoint_value
    )
    scales = tuple(right.scale(degree) for degree in range(degree_max + 1))

    rows: list[BorelPrefixRightJetDegreeCertificate] = []
    for degree in range(degree_max + 1):
        prefix_scales = scales[: degree + 1]
        radius = Fraction(1, 2 * max(prefix_scales))
        row = BorelPrefixRightJetDegreeCertificate(
            degree=degree,
            point=point,
            endpoint=endpoint_value,
            prefix_scales=prefix_scales,
            common_plateau_radius=radius,
            prescribed_normal_jet=normal_jets[degree],
            prefix_endpoint_time_jet=_prefix_endpoint_derivative(normal_jets, degree),
        )
        if not row.finite_prefix_right_jet_verified:
            raise ArithmeticError("finite Borel prefix right-jet invariant failed")
        rows.append(row)

    result = Section10BorelPrefixRightJetsCertificate(
        point=point,
        endpoint=endpoint_value,
        max_degree=degree_max,
        rows=tuple(rows),
    )
    if not result.formal_prefix_ready:
        raise ArithmeticError("Section 10 Borel prefix right-jet certificate invariant failed")
    return result
