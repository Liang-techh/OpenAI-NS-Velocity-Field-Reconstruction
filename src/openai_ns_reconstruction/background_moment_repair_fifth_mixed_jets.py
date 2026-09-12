"""Analytic fifth-mixed U-jet bridge for the Lemma 5.2 compact repair.

The landed Eq. (5.2) fourth-mixed regular-flux adapter needs exactly three
additional total-order-five axial derivatives:
``U_XXetaetaeta``, ``U_Xetaetaetaeta`` and ``U_etaetaetaetaeta``.
The compact Lemma 5.2 repair already has arbitrary finite eta jets of its
coefficient functions, while its radial bumps already expose exact X jets
through second order.

This module composes those pieces without finite differencing in production.
The unrepaired fifth-mixed U jet and the moment/patch eta jets remain explicit
upstream inputs, so the result is fail-closed formal-structure infrastructure,
not a paper-exact recursive coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_fourth_mixed_jets import (
    Lemma52RepairedFourthMixedJetAdapter,
)
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_second_jets import _bump_x_second_jet
from .background_omega_second_parameter_jet import RegularFluxFourthMixedJet
from .background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
    regular_flux_fourth_mixed_jet_eq_5_2,
)

BaseFifthMixedJetProvider = Callable[[float, float], AxialFifthMixedJet]


def _provider_jet(
    provider: BaseFifthMixedJetProvider,
    X: float,
    eta: float,
) -> AxialFifthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialFifthMixedJet):
        raise TypeError("base_U_fifth_mixed_jet must return AxialFifthMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedFifthMixedJetAdapter:
    """Lift one compact Lemma 5.2 U repair to the fifth-mixed Eq. (5.2) jet.

    For ``Delta U=sum_j alpha_j(eta)b_j(sqrt(2X))``, the only new fields beyond
    the landed repaired fourth-mixed adapter are

    ``U_XXetaetaeta    += sum_j alpha_j'''   d_X^2 b_j``,
    ``U_Xetaetaetaeta  += sum_j alpha_j''''  d_X b_j``, and
    ``U_etaetaetaetaeta += sum_j alpha_j''''' b_j``.

    The first twelve fields are obtained from the existing fourth-mixed repair
    using the exact projection of the same base fifth-mixed provider.  Thus the
    stronger jet cannot silently diverge from the already-landed lower layer.
    """

    profile: Lemma52RepairedProfileAdapter
    base_U_fifth_mixed_jet: BaseFifthMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired axial fifth-mixed jet from caller-supplied unrepaired "
        "fifth-mixed U data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        if not callable(self.base_U_fifth_mixed_jet):
            raise TypeError("base_U_fifth_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    def _fourth_bridge(self) -> Lemma52RepairedFourthMixedJetAdapter:
        return Lemma52RepairedFourthMixedJetAdapter(
            self.profile,
            lambda X, eta: _provider_jet(
                self.base_U_fifth_mixed_jet, X, eta
            ).fourth(),
        )

    def U_fifth_mixed_jet(self, X: float, eta: float) -> AxialFifthMixedJet:
        """Return the compact-repaired axial fifth-mixed jet."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_U_fifth_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.u_bumps)
        if not self._has_active_bump(X, bumps):
            # The compact C-infinity repair is exactly flat off support, so no
            # fifth eta moment/patch row is needed there.
            return base

        lower = self._fourth_bridge().U_fourth_mixed_jet(X, eta)
        coeff = self.profile.repair_coefficients(eta, max_order=5)
        rows = coeff.alpha_derivatives
        radial_jets = tuple(_bump_x_second_jet(bump, X) for bump in bumps)

        def correction(eta_order: int, radial_order: int) -> float:
            return math.fsum(
                float(rows[eta_order, j]) * radial_jets[j][radial_order]
                for j in range(len(bumps))
            )

        return AxialFifthMixedJet(
            value=lower.value,
            radial=lower.radial,
            radial2=lower.radial2,
            parameter=lower.parameter,
            radial_parameter=lower.radial_parameter,
            parameter2=lower.parameter2,
            radial2_parameter=lower.radial2_parameter,
            radial_parameter2=lower.radial_parameter2,
            parameter3=lower.parameter3,
            radial2_parameter2=lower.radial2_parameter2,
            radial_parameter3=lower.radial_parameter3,
            parameter4=lower.parameter4,
            radial2_parameter3=base.radial2_parameter3 + correction(3, 2),
            radial_parameter4=base.radial_parameter4 + correction(4, 1),
            parameter5=base.parameter5 + correction(5, 0),
        )

    def beta_fourth_mixed_jet(
        self,
        h: float,
        order: int,
        X: float,
        eta: float,
        *,
        quadrature_points: int = 32,
    ) -> RegularFluxFourthMixedJet:
        """Feed the repaired fifth-mixed U jet into the analytic Eq. (5.2) map."""

        return regular_flux_fourth_mixed_jet_eq_5_2(
            h,
            order,
            X,
            eta,
            self.U_fifth_mixed_jet,
            quadrature_points=quadrature_points,
        )
