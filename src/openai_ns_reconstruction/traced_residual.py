"""Value-level endpoint trace and Borel glue for Section 10.

The pinned OpenAI formalization defines ``CandidateFromLimits.tracedResidual``
as ``SpacetimeEndpoint.extendTrace`` applied to the actual closed-past
Navier--Stokes residual: for ``t < T`` it is the past residual, while for
``t >= T`` its value is the degree-zero endpoint tensor.  The final force then
uses ``SpacetimeGluing.smoothExtension`` to keep that closed-past trace and
realize the normal endpoint jets on the future side with a Taylor--Borel sum.

This module makes only that *value-level algebra* executable.  The supplied
full spacetime jets are not evidence that they are locally-uniform limits of
actual residual derivatives, and a supplied scale is not an analytic
majorant certificate.  Consequently these helpers remain formal-structure
infrastructure and must not be used to promote the paper-exact gate.
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

from .coordinates import _finite
from .endpoint_borel import LocalScale
from .endpoint_jets import (
    FullSpacetimeJetFamily,
    borel_extension_from_full_spacetime_jets,
    validate_full_spacetime_jet,
)
from .verify import finite_vector

VectorField = Callable[[float, float, float, float], np.ndarray]


def traced_residual_from_full_spacetime_jets(
    past_residual: VectorField,
    full_jets: FullSpacetimeJetFamily,
    *,
    endpoint: float = 1.0,
) -> VectorField:
    """Executable counterpart of ``CandidateFromLimits.tracedResidual``.

    The branch is exactly the pinned ``SpacetimeEndpoint.extendTrace`` value
    rule.  The past branch is short-circuited, so evaluating ``t >= endpoint``
    does not sample an unresolved past field.  Conversely, evaluating
    ``t < endpoint`` does not demand any endpoint-jet data.

    Only the order-zero tensor is used on/after the endpoint because this
    helper constructs the *trace value*.  Higher tensors enter through the
    Taylor--Borel future extension in
    :func:`formal_force_from_full_spacetime_jets`.
    """
    if not callable(past_residual):
        raise ValueError("past_residual must be callable")
    if not callable(full_jets):
        raise ValueError("full_jets must be callable")
    endpoint = _finite(endpoint, "endpoint")

    def traced(x: float, y: float, z: float, t: float) -> np.ndarray:
        x = _finite(x, "x")
        y = _finite(y, "y")
        z = _finite(z, "z")
        t = _finite(t, "t")
        if t < endpoint:
            return finite_vector(past_residual(x, y, z, t))
        tensor0 = validate_full_spacetime_jet(full_jets(0, x, y, z), 0)
        return finite_vector(tensor0)

    return traced


def formal_force_from_full_spacetime_jets(
    past_residual: VectorField,
    full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    *,
    endpoint: float = 1.0,
) -> VectorField:
    """Compose the endpoint trace with the landed normal-jet Borel branch.

    This mirrors the *branch structure* of ``CandidateFromLimits.force``:

    * ``t <= T`` uses the traced closed-past residual (with the endpoint value
      filled by the order-zero limit tensor), and
    * ``t > T`` uses the Taylor--Borel right extension of the repeated-time
      contractions of the full endpoint tensors.

    The returned field is deliberately named ``formal_force`` rather than
    ``force``: no smoothness, derivative matching, locally-uniform convergence,
    or paper-exact claim follows merely from supplying arrays and scales.
    """
    traced = traced_residual_from_full_spacetime_jets(
        past_residual, full_jets, endpoint=endpoint
    )
    right = borel_extension_from_full_spacetime_jets(
        full_jets, local_scale, endpoint=endpoint
    )
    return right.glue_to_past(traced)
