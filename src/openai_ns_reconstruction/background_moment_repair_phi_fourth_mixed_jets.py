"""Analytic fourth-mixed ``phi_n`` jets for the Lemma 5.2 compact swirl repair.

The hierarchy-owned second-eta strict-lower source needs the second eta jet of
its angular preceding-diffusion term.  The landed analytic operator for that
term consumes ``ProfileFourthMixedJet``, while the compact Lemma 5.2 phi bridge
previously stopped at third-mixed data.

On each active E-repair support,

    Delta phi_n = C sum_j beta_j(eta) b_j(R) / R,   R = sqrt(2 X).

The radial factors are eta-independent.  Thus the three new total-order-four
entries pair beta'' with ``d_X^2(b/R)``, beta''' with ``d_X(b/R)``, and beta''''
with ``b/R``.  Off compact support the supplied base fourth-mixed jet is
returned exactly, and no fourth eta row is requested.

The unrepaired fourth-mixed phi jet, normalization C, and moment/patch eta jets
remain upstream inputs.  This is Stage-2 ``formal-structure`` infrastructure;
it does not make caller data or finite-order test fixtures paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_phi_jets import _bump_over_radius_x_second_jet
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
)


BasePhiFourthMixedJetProvider = Callable[[float, float], ProfileFourthMixedJet]


def _provider_jet(
    provider: BasePhiFourthMixedJetProvider,
    X: float,
    eta: float,
) -> ProfileFourthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, ProfileFourthMixedJet):
        raise TypeError("base_phi_fourth_mixed_jet must return ProfileFourthMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedPhiFourthMixedJetAdapter:
    """Lift the compact repaired ``phi_n`` profile to a fourth-mixed scalar jet."""

    profile: Lemma52RepairedProfileAdapter
    C: float
    base_phi_fourth_mixed_jet: BasePhiFourthMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired phi=C E/R fourth-mixed bridge from caller-supplied "
        "unrepaired phi data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        C = float(self.C)
        if not math.isfinite(C) or C <= 0.0:
            raise ValueError("C must be finite and positive")
        if not callable(self.base_phi_fourth_mixed_jet):
            raise TypeError("base_phi_fourth_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")
        object.__setattr__(self, "C", C)

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    def phi_fourth_mixed_jet(self, X: float, eta: float) -> ProfileFourthMixedJet:
        """Return the repaired PositiveAxis ``phi_n`` fourth-mixed jet."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_phi_fourth_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.e_bumps)
        if not self._has_active_bump(X, bumps):
            return base

        coeff = self.profile.repair_coefficients(eta, max_order=4)
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
                    "compact phi repair fourth-mixed correction is outside binary64 range"
                )
            return value

        return ProfileFourthMixedJet(
            value=base.value + correction(0, 0),
            radial=base.radial + correction(0, 1),
            radial2=base.radial2 + correction(0, 2),
            parameter=base.parameter + correction(1, 0),
            radial_parameter=base.radial_parameter + correction(1, 1),
            parameter2=base.parameter2 + correction(2, 0),
            radial2_parameter=base.radial2_parameter + correction(1, 2),
            radial_parameter2=base.radial_parameter2 + correction(2, 1),
            parameter3=base.parameter3 + correction(3, 0),
            radial2_parameter2=base.radial2_parameter2 + correction(2, 2),
            radial_parameter3=base.radial_parameter3 + correction(3, 1),
            parameter4=base.parameter4 + correction(4, 0),
        )
