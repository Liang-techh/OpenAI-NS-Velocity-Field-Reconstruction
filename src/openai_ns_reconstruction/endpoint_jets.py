"""Typed endpoint-jet adapters for the Section 10 Taylor--Borel glue.

The pinned Lean construction exposes two related interfaces:

* ``CandidateFromLimits`` assumes locally uniform limits ``L x n`` of the
  *full* order-n spacetime Frechet derivatives of the actual Navier--Stokes
  residual;
* ``SpacetimeGluing`` feeds repeated derivatives in its fixed time direction
  ``timeVector = (1, 0)`` to ``SpatialBorelExtension.rightExtension``.

For an executable dense tensor representation, this module performs exactly
that contraction.  A degree-n tensor has shape ``(3,) + (4,)*n``: the leading
axis is the vector-field output and each derivative slot is ordered
``(time, x, y, z)``.  The normal time jet is therefore ``J[:,0,...,0]``.

This is an interface adapter, not an endpoint regularity certificate.  Caller
supplied tensors are not evidence that they are limits of derivatives of the
actual residual, that the convergence is locally uniform, or that the
all-order derivative-majorant scale hypotheses used by the Lean Borel proof
hold.  Consequently this module remains formal-structure infrastructure and
must not be used to promote the paper-exact gate.
"""
from __future__ import annotations

from collections.abc import Callable
import operator

import numpy as np

from .endpoint_borel import BorelRightExtension, LocalScale
from .verify import finite_vector

FullSpacetimeJetFamily = Callable[[int, float, float, float], object]
SpatialJetFamily = Callable[[int, float, float, float], np.ndarray]


def _degree(value: object) -> int:
    """Validate an exact nonnegative integer derivative degree."""
    if isinstance(value, bool):
        raise ValueError("degree must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError("degree must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError("degree must be a nonnegative integer")
    return int(out)


def validate_full_spacetime_jet(value: object, degree: int) -> np.ndarray:
    """Validate the dense representation of an order-``degree`` full jet.

    The output shape is ``(3,) + (4,)*degree``.  This convention represents a
    vector-valued multilinear map on spacetime with ordered coordinate basis
    ``(time, x, y, z)``.  No symmetry between derivative slots is assumed.
    """
    degree = _degree(degree)
    tensor = np.asarray(value, dtype=float)
    expected = (3,) + (4,) * degree
    if tensor.shape != expected:
        raise ValueError(
            f"order-{degree} full spacetime jet must have shape {expected}, got {tensor.shape}"
        )
    if not np.all(np.isfinite(tensor)):
        raise ValueError("full spacetime jet must contain only finite values")
    return tensor


def normal_time_jet_from_full(
    full_jets: FullSpacetimeJetFamily,
    degree: int,
    x: float,
    y: float,
    z: float,
) -> np.ndarray:
    """Contract every derivative slot of the full jet with the time basis.

    In the pinned ``SpacetimeGluing.lean`` this is repeated differentiation in
    ``timeVector = (1,0)``.  In the dense coordinate convention used here, the
    contraction is exactly the all-zero derivative multi-index.
    """
    if not callable(full_jets):
        raise ValueError("full_jets must be callable")
    degree = _degree(degree)
    tensor = validate_full_spacetime_jet(full_jets(degree, x, y, z), degree)
    index = (slice(None),) + (0,) * degree
    return finite_vector(tensor[index])


def normal_jet_family(full_jets: FullSpacetimeJetFamily) -> SpatialJetFamily:
    """Adapt full endpoint derivative tensors to the Borel normal-jet API."""
    if not callable(full_jets):
        raise ValueError("full_jets must be callable")

    def normal(degree: int, x: float, y: float, z: float) -> np.ndarray:
        return normal_time_jet_from_full(full_jets, degree, x, y, z)

    return normal


def borel_extension_from_full_spacetime_jets(
    full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    *,
    endpoint: float = 1.0,
) -> BorelRightExtension:
    """Compose the full-jet contraction with the existing Borel evaluator.

    This helper deliberately accepts the same caller-supplied ``local_scale``
    boundary as ``BorelRightExtension``.  It does not construct the official
    derivative-majorant scales and does not assert that ``full_jets`` arise
    from actual residual limits.
    """
    return BorelRightExtension(
        normal_jet_family(full_jets), local_scale, endpoint=endpoint
    )
