"""Analytic fifth-mixed ``phi_n`` jets for the Lemma 5.2 compact swirl repair.

The hierarchy-owned third-eta strict-lower source will need one derivative tier
above the landed fourth-mixed ``phi_n`` bridge.  On each active E-repair support,

    Delta phi_n = C sum_j beta_j(eta) b_j(R) / R,   R = sqrt(2 X).

The radial factors are eta-independent.  Therefore the only entries beyond the
landed fourth-mixed bridge are ``phi_XXetaetaeta``, ``phi_Xetaetaetaeta`` and
``phi_etaetaetaetaeta``.  They use beta''', beta'''' and beta''''' together
with the already-landed exact X-second jet of ``b/R``.

Off compact support the supplied base fifth-mixed jet is returned exactly and
no fifth eta row is requested.  The unrepaired fifth-mixed phi jet,
normalization C, and moment/patch eta jets remain upstream inputs.  This is
Stage-2 ``formal-structure`` infrastructure, not a paper-exact coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_phi_fourth_mixed_jets import (
    Lemma52RepairedPhiFourthMixedJetAdapter,
)
from .background_moment_repair_phi_jets import _bump_over_radius_x_second_jet
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
)


@dataclass(frozen=True)
class ProfileFifthMixedJet:
    """Scalar X/eta jet through the three mixed fields of total order five."""

    value: float
    radial: float
    radial2: float
    parameter: float
    radial_parameter: float
    parameter2: float
    radial2_parameter: float
    radial_parameter2: float
    parameter3: float
    radial2_parameter2: float
    radial_parameter3: float
    parameter4: float
    radial2_parameter3: float
    radial_parameter4: float
    parameter5: float

    def fourth(self) -> ProfileFourthMixedJet:
        """Project exactly to the landed total-order-four scalar jet."""

        return ProfileFourthMixedJet(
            value=self.value,
            radial=self.radial,
            radial2=self.radial2,
            parameter=self.parameter,
            radial_parameter=self.radial_parameter,
            parameter2=self.parameter2,
            radial2_parameter=self.radial2_parameter,
            radial_parameter2=self.radial_parameter2,
            parameter3=self.parameter3,
            radial2_parameter2=self.radial2_parameter2,
            radial_parameter3=self.radial_parameter3,
            parameter4=self.parameter4,
        )


BasePhiFifthMixedJetProvider = Callable[[float, float], ProfileFifthMixedJet]


def _provider_jet(
    provider: BasePhiFifthMixedJetProvider,
    X: float,
    eta: float,
) -> ProfileFifthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, ProfileFifthMixedJet):
        raise TypeError("base_phi_fifth_mixed_jet must return ProfileFifthMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedPhiFifthMixedJetAdapter:
    """Lift the compact repaired ``phi_n`` profile to a fifth-mixed scalar jet.

    The first twelve fields are obtained from the existing fourth-mixed repair
    using the exact projection of the same base fifth-mixed provider.  Thus the
    stronger jet cannot silently diverge from the already-landed lower layer.
    """

    profile: Lemma52RepairedProfileAdapter
    C: float
    base_phi_fifth_mixed_jet: BasePhiFifthMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired phi=C E/R fifth-mixed bridge from caller-supplied "
        "unrepaired phi data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        C = float(self.C)
        if not math.isfinite(C) or C <= 0.0:
            raise ValueError("C must be finite and positive")
        if not callable(self.base_phi_fifth_mixed_jet):
            raise TypeError("base_phi_fifth_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")
        object.__setattr__(self, "C", C)

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    def _fourth_bridge(self) -> Lemma52RepairedPhiFourthMixedJetAdapter:
        return Lemma52RepairedPhiFourthMixedJetAdapter(
            self.profile,
            self.C,
            lambda X, eta: _provider_jet(
                self.base_phi_fifth_mixed_jet, X, eta
            ).fourth(),
        )

    def phi_fifth_mixed_jet(self, X: float, eta: float) -> ProfileFifthMixedJet:
        """Return the repaired PositiveAxis ``phi_n`` fifth-mixed jet."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_phi_fifth_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.e_bumps)
        if not self._has_active_bump(X, bumps):
            return base

        lower = self._fourth_bridge().phi_fourth_mixed_jet(X, eta)
        coeff = self.profile.repair_coefficients(eta, max_order=5)
        rows = coeff.beta_derivatives
        radial_jets = tuple(
            _bump_over_radius_x_second_jet(bump, X) for bump in bumps
        )

        def correction(eta_order: int, radial_order: int) -> float:
            value = self.C * math.fsum(
                float(rows[eta_order, j]) * radial_jets[j][radial_order]
                for j in range(len(bumps))
            )
            if not math.isfinite(value):
                raise OverflowError(
                    "compact phi repair fifth-mixed correction is outside binary64 range"
                )
            return value

        return ProfileFifthMixedJet(
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
