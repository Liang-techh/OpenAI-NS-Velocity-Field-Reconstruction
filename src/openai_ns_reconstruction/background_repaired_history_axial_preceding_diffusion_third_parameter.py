"""Hierarchy-owned axial third-eta preceding diffusion for Section 5.

PR #280 landed the exact analytic third-eta preceding-diffusion primitive and a
hierarchy-owned angular bridge.  The strict-lower axial row uses the same
operator on ``U_(n-1)`` with paper power ``-1/2-h``.  The repaired lower-history
hierarchy already owns the required fifth-mixed ``U`` jet, so this module closes
that remaining preceding-diffusion ownership seam without a caller-supplied
profile/source jet.

Production performs no finite differencing, fitting, generic-cutoff
substitution, or sampled derivative reconstruction.  Genuine leading strong
data remain an Issue-#1 dependency, so this is Stage-2 ``formal-structure``
solver infrastructure rather than a paper-exact coefficient or convergence
claim.
"""
from __future__ import annotations

from numbers import Integral

from .background_moment_repair_phi_fifth_mixed_jets import ProfileFifthMixedJet
from .background_preceding_diffusion_third_parameter_jet import (
    PrecedingDiffusionThirdParameterJet,
    preceding_diffusion_third_parameter_jet,
)
from .background_regular_flux_fourth_mixed_jets import AxialFifthMixedJet
from .background_repaired_history_fifth_mixed import (
    Section5LowerHistoryFifthMixedHierarchy,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _profile_fifth_from_axial(jet: AxialFifthMixedJet) -> ProfileFifthMixedJet:
    """Re-type the common mixed-jet rows without changing any values."""

    if not isinstance(jet, AxialFifthMixedJet):
        raise TypeError("jet must be an AxialFifthMixedJet")
    return ProfileFifthMixedJet(
        value=jet.value,
        radial=jet.radial,
        radial2=jet.radial2,
        parameter=jet.parameter,
        radial_parameter=jet.radial_parameter,
        parameter2=jet.parameter2,
        radial2_parameter=jet.radial2_parameter,
        radial_parameter2=jet.radial_parameter2,
        parameter3=jet.parameter3,
        radial2_parameter2=jet.radial2_parameter2,
        radial_parameter3=jet.radial_parameter3,
        parameter4=jet.parameter4,
        radial2_parameter3=jet.radial2_parameter3,
        radial_parameter4=jet.radial_parameter4,
        parameter5=jet.parameter5,
    )


def hierarchy_owned_axial_preceding_diffusion_third_parameter_jet(
    hierarchy: Section5LowerHistoryFifthMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> PrecedingDiffusionThirdParameterJet:
    """Return the hierarchy-owned axial preceding-diffusion third eta jet.

    ``order`` is the positive coefficient order being solved.  The paper's
    strict-lower axial preceding-diffusion row uses ``U_(order-1)`` with axial
    power ``-1/2-h``.  The authoritative repaired fifth-mixed U provider is
    queried before the operator is evaluated, so missing strong history fails
    closed even at the axis.
    """

    if not isinstance(hierarchy, Section5LowerHistoryFifthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryFifthMixedHierarchy"
        )
    order = _positive_order(order)
    if len(hierarchy.sources) < order:
        raise ValueError(
            "requested source order is missing one or more strict-lower coefficients"
        )

    axial = hierarchy.axial_fifth_mixed_jet(order - 1, X, eta)
    return preceding_diffusion_third_parameter_jet(
        hierarchy.h,
        -0.5 - hierarchy.h,
        order,
        X,
        eta,
        _profile_fifth_from_axial(axial),
    )
