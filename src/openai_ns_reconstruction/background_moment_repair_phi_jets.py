"""Analytic ``phi_n`` second jets for the Lemma 5.2 compact swirl repair.

Lemma 5.2 repairs the angular coefficient ``E_n`` by compactly supported bumps
and then defines the PositiveAxis unknown by the exact paper identity

    phi_n = C E_n / R,    R = sqrt(2 X).

The landed second-jet repair already differentiates the compact ``E_n`` bumps,
but the strict-lower-history solver consumes second ``(X, eta)`` jets of
``phi_n`` rather than ``E_n``.  This module applies the ``C/R`` factor
analytically on the positive-radius repair supports.  Off those supports the
repair is identically flat, so the supplied inner/base ``phi_n`` jet is returned
unchanged, including at the axis; production never divides sampled axis data by
``R`` and never finite-differences derivatives.

The base ``phi_n`` jet, the manuscript normalization ``C``, and the moment/patch
eta-jets remain upstream inputs.  Consequently this is Section-5 solver
infrastructure / ``formal-structure`` only and must not be promoted to a
paper-exact recursive coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_lower_history_source import ProfileSecondJet
from .background_moment_repair import CompactMomentBump
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_second_jets import _bump_x_second_jet


BasePhiSecondJetProvider = Callable[[float, float], ProfileSecondJet]


def _bump_over_radius_x_second_jet(
    bump: CompactMomentBump,
    X: float,
) -> tuple[float, float, float]:
    """Return ``(g, d_X g, d_X^2 g)`` for ``g=b(sqrt(2X))/sqrt(2X)``.

    Every Lemma-5.2 repair bump is supported on a strictly positive radial
    interval.  Therefore division by ``R`` occurs only where ``R>0``.  Outside
    (and at the flat support boundary) all three derivatives are exactly zero.
    """

    X = float(X)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not isinstance(bump, CompactMomentBump):
        raise TypeError("bump must be a CompactMomentBump")

    radius = math.sqrt(2.0 * X)
    if radius <= bump.left or radius >= bump.right:
        return 0.0, 0.0, 0.0

    value, d_X, d_X2 = _bump_x_second_jet(bump, X)
    inv_r = 1.0 / radius
    g = value * inv_r
    g_X = d_X * inv_r - value * inv_r**3
    g_X2 = d_X2 * inv_r - 2.0 * d_X * inv_r**3 + 3.0 * value * inv_r**5
    if not all(math.isfinite(v) for v in (g, g_X, g_X2)):
        raise OverflowError("compact phi repair second jet is outside binary64 range")
    return g, g_X, g_X2


def _provider_jet(
    provider: BasePhiSecondJetProvider,
    X: float,
    eta: float,
) -> ProfileSecondJet:
    value = provider(X, eta)
    if not isinstance(value, ProfileSecondJet):
        raise TypeError("base_phi_second_jet must return ProfileSecondJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedPhiSecondJetAdapter:
    """Lift the compact ``E_n`` repair to the PositiveAxis ``phi_n`` jet.

    On the repair support,

    ``Delta phi_n = C sum_j beta_j(eta) b_j(R) / R``.

    The eta derivatives act only on ``beta_j`` and the X derivatives act only
    on ``b_j(R)/R``.  The returned :class:`ProfileSecondJet` is therefore the
    exact compact-repair contribution added to the supplied inner/base
    ``phi_n`` jet.
    """

    profile: Lemma52RepairedProfileAdapter
    C: float
    base_phi_second_jet: BasePhiSecondJetProvider
    provenance: str = (
        "Lemma 5.2 repaired phi=C E/R second-jet bridge from caller-supplied "
        "inner phi data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        C = float(self.C)
        if not math.isfinite(C) or C <= 0.0:
            raise ValueError("C must be finite and positive")
        if not callable(self.base_phi_second_jet):
            raise TypeError("base_phi_second_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")
        object.__setattr__(self, "C", C)

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    def phi_second_jet(self, X: float, eta: float) -> ProfileSecondJet:
        """Return the repaired PositiveAxis ``phi_n`` second jet at ``(X,eta)``."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_phi_second_jet, X, eta)
        bumps = tuple(self.profile.repair.e_bumps)
        if not self._has_active_bump(X, bumps):
            # The compact correction is identically flat here.  In particular,
            # at X=0 no C/R division is ever evaluated.
            return base

        coeff = self.profile.repair_coefficients(eta, max_order=2)
        rows = coeff.beta_derivatives
        radial_jets = tuple(_bump_over_radius_x_second_jet(bump, X) for bump in bumps)

        def correction(eta_order: int, radial_order: int) -> float:
            return self.C * math.fsum(
                float(rows[eta_order, j]) * radial_jets[j][radial_order]
                for j in range(len(bumps))
            )

        return ProfileSecondJet(
            value=base.value + correction(0, 0),
            radial=base.radial + correction(0, 1),
            radial2=base.radial2 + correction(0, 2),
            parameter=base.parameter + correction(1, 0),
            radial_parameter=base.radial_parameter + correction(1, 1),
            parameter2=base.parameter2 + correction(2, 0),
        )
