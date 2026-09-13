"""Analytic third eta jet for the Section 5 preceding-diffusion term.

The strict-lower PositiveAxis source contains

    precedingDiffusion = Z_(b-D) (Z_b F_(n-1)).

The landed second-parameter bridge owns the value, first eta derivative and
second eta derivative.  Differentiating once more requires exactly the three
eta-bearing total-order-five entries carried by ``ProfileFifthMixedJet``:
``F_XXetaetaeta``, ``F_Xetaetaetaeta`` and ``F_etaetaetaetaeta``.

This module computes only the new third derivative analytically and delegates
all lower derivatives to the landed second-parameter implementation.  It also
provides a narrow hierarchy-owned angular bridge: for positive coefficient
order ``n`` it consumes the strong repaired hierarchy's authoritative
``phi_(n-1)`` fifth-mixed jet.  No production finite difference, generic cutoff,
sampled derivative table, or fitted coefficient is admitted.

The hierarchy still depends on genuine Issue-#1 leading strong data and on the
explicit upstream contracts used by the Lemma 5.2 repair.  This is therefore
Stage-2 ``formal-structure`` solver infrastructure, not a paper-exact velocity
or a convergence claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

from .background_moment_repair_phi_fifth_mixed_jets import ProfileFifthMixedJet
from .background_preceding_diffusion_second_parameter_jet import (
    preceding_diffusion_second_parameter_jet,
)
from .background_repaired_history_phi_fifth_mixed import (
    Section5LowerHistoryPhiFifthMixedHierarchy,
)
from .coordinates import validate_h


@dataclass(frozen=True)
class PrecedingDiffusionThirdParameterJet:
    """Value and the first three analytic eta derivatives."""

    value: float
    parameter: float
    parameter2: float
    parameter3: float

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2", "parameter3"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _new_fifth_entries(jet: ProfileFifthMixedJet) -> None:
    for name in ("radial2_parameter3", "radial_parameter4", "parameter5"):
        value = float(getattr(jet, name))
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite")


def _inverse_ell_jet(h: float, eta: float) -> tuple[float, float, float, float, float]:
    """Return eta derivatives 0..4 of ``1/(1-2 h eta^2)``."""

    ell = 1.0 - 2.0 * h * eta * eta
    s0 = 1.0 / ell
    s1 = 4.0 * h * eta / ell**2
    s2 = 4.0 * h / ell**2 + 32.0 * h * h * eta * eta / ell**3
    s3 = 96.0 * h * h * eta / ell**3 + 384.0 * h**3 * eta**3 / ell**4
    s4 = (
        96.0 * h * h / ell**3
        + 2304.0 * h**3 * eta * eta / ell**4
        + 6144.0 * h**4 * eta**4 / ell**5
    )
    return s0, s1, s2, s3, s4


def _quotient_derivatives(
    numerator: tuple[float, ...],
    inverse_ell: tuple[float, ...],
    max_order: int,
) -> tuple[float, ...]:
    """Leibniz jets for ``numerator / ell`` from exact derivative rows."""

    return tuple(
        math.fsum(
            math.comb(order, derivative)
            * numerator[derivative]
            * inverse_ell[order - derivative]
            for derivative in range(order + 1)
        )
        for order in range(max_order + 1)
    )


def preceding_diffusion_third_parameter_jet(
    h: float,
    power: float,
    order: int,
    X: float,
    eta: float,
    jet: ProfileFifthMixedJet,
) -> PrecedingDiffusionThirdParameterJet:
    """Return the first three analytic eta derivatives of preceding diffusion.

    The value/first/second rows are delegated to the landed
    :func:`preceding_diffusion_second_parameter_jet` after exact projection of
    the stronger fifth-mixed input.  Only ``partial_eta^3`` is assembled here.
    """

    if not isinstance(jet, ProfileFifthMixedJet):
        raise TypeError("jet must be ProfileFifthMixedJet")
    lower = preceding_diffusion_second_parameter_jet(
        h,
        power,
        order,
        X,
        eta,
        jet.fourth(),
    )
    _new_fifth_entries(jet)

    h = validate_h(h)
    order = _positive_order(order)
    power = float(power)
    X = float(X)
    eta = float(eta)
    if not math.isfinite(power):
        raise ValueError("power must be finite")
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")

    D = 0.5 - h
    d = 1.0 - eta * eta
    b = power + 2.0 * (order - 1) * h

    # N is the numerator of the inner Z_b F.  N^(0)..N^(4) and
    # N_X^(0)..N_X^(3) are written explicitly so the three new fifth-mixed
    # entries enter only where the analytic differentiation requires them.
    N0 = 2.0 * eta * b * jet.value + d * jet.parameter - 2.0 * eta * X * jet.radial
    N1 = (
        2.0 * b * jet.value
        + 2.0 * eta * (b - 1.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    N2 = (
        (4.0 * b - 2.0) * jet.parameter
        + 2.0 * eta * (b - 2.0) * jet.parameter2
        + d * jet.parameter3
        - 4.0 * X * jet.radial_parameter
        - 2.0 * eta * X * jet.radial_parameter2
    )
    N3 = (
        6.0 * (b - 1.0) * jet.parameter2
        + 2.0 * eta * (b - 3.0) * jet.parameter3
        + d * jet.parameter4
        - 6.0 * X * jet.radial_parameter2
        - 2.0 * eta * X * jet.radial_parameter3
    )
    N4 = (
        (8.0 * b - 12.0) * jet.parameter3
        + 2.0 * eta * (b - 4.0) * jet.parameter4
        + d * jet.parameter5
        - 8.0 * X * jet.radial_parameter3
        - 2.0 * eta * X * jet.radial_parameter4
    )

    NX0 = (
        2.0 * eta * (b - 1.0) * jet.radial
        + d * jet.radial_parameter
        - 2.0 * eta * X * jet.radial2
    )
    NX1 = (
        2.0 * (b - 1.0) * jet.radial
        + 2.0 * eta * (b - 2.0) * jet.radial_parameter
        + d * jet.radial_parameter2
        - 2.0 * X * jet.radial2
        - 2.0 * eta * X * jet.radial2_parameter
    )
    NX2 = (
        (4.0 * b - 6.0) * jet.radial_parameter
        + 2.0 * eta * (b - 3.0) * jet.radial_parameter2
        + d * jet.radial_parameter3
        - 4.0 * X * jet.radial2_parameter
        - 2.0 * eta * X * jet.radial2_parameter2
    )
    NX3 = (
        6.0 * (b - 2.0) * jet.radial_parameter2
        + 2.0 * eta * (b - 4.0) * jet.radial_parameter3
        + d * jet.radial_parameter4
        - 6.0 * X * jet.radial2_parameter2
        - 2.0 * eta * X * jet.radial2_parameter3
    )

    inverse_ell = _inverse_ell_jet(h, eta)
    first = _quotient_derivatives((N0, N1, N2, N3, N4), inverse_ell, 4)
    first_X = _quotient_derivatives((NX0, NX1, NX2, NX3), inverse_ell, 3)

    c = b - D
    Q0 = 2.0 * eta * c * first[0] + d * first[1] - 2.0 * eta * X * first_X[0]
    Q1 = (
        2.0 * c * first[0]
        + 2.0 * eta * (c - 1.0) * first[1]
        + d * first[2]
        - 2.0 * X * first_X[0]
        - 2.0 * eta * X * first_X[1]
    )
    Q2 = (
        (4.0 * c - 2.0) * first[1]
        + 2.0 * eta * (c - 2.0) * first[2]
        + d * first[3]
        - 4.0 * X * first_X[1]
        - 2.0 * eta * X * first_X[2]
    )
    Q3 = (
        6.0 * (c - 1.0) * first[2]
        + 2.0 * eta * (c - 3.0) * first[3]
        + d * first[4]
        - 6.0 * X * first_X[2]
        - 2.0 * eta * X * first_X[3]
    )

    s0, s1, s2, s3, _s4 = inverse_ell
    parameter3 = Q3 * s0 + 3.0 * Q2 * s1 + 3.0 * Q1 * s2 + Q0 * s3
    if not math.isfinite(parameter3):
        raise OverflowError(
            "preceding diffusion third parameter derivative is outside binary64 range"
        )

    return PrecedingDiffusionThirdParameterJet(
        value=lower.value,
        parameter=lower.parameter,
        parameter2=lower.parameter2,
        parameter3=parameter3,
    )


def hierarchy_owned_angular_preceding_diffusion_third_parameter_jet(
    hierarchy: Section5LowerHistoryPhiFifthMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> PrecedingDiffusionThirdParameterJet:
    """Return the hierarchy-owned angular preceding-diffusion third eta jet.

    ``order`` is the positive coefficient order being solved.  The paper's
    strict-lower preceding-diffusion row uses ``phi_(order-1)`` with angular
    power ``-1-h``.  The fifth-mixed hierarchy query is deliberately performed
    before the operator evaluation so missing genuine strong data fail closed,
    including on the axis.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiFifthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiFifthMixedHierarchy"
        )
    order = _positive_order(order)
    if len(hierarchy.sources) < order:
        raise ValueError(
            "requested source order is missing one or more strict-lower coefficients"
        )
    jet = hierarchy.phi_fifth_mixed_jet(order - 1, X, eta)
    return preceding_diffusion_third_parameter_jet(
        hierarchy.h,
        -1.0 - hierarchy.h,
        order,
        X,
        eta,
        jet,
    )
