"""Analytic third-mixed ``phi_n`` jets for the Lemma 5.2 compact swirl repair.

The strict-lower PositiveAxis source now has an analytic eta derivative for the
``precedingDiffusion = Z_(b-D)(Z_b F_(n-1))`` term. That derivative consumes a
scalar third-mixed jet, while the landed Lemma 5.2 ``phi_n`` bridge only owns a
second jet. This module lifts the same exact compact repair one derivative level
higher.

On every active E-repair support,

    Delta phi_n = C sum_j beta_j(eta) b_j(R) / R,   R = sqrt(2 X).

Thus the three additional derivatives pair beta'_j with ``d_X^2(b_j/R)``,
beta''_j with ``d_X(b_j/R)``, and beta'''_j with ``b_j/R``. The repair supports
are strictly away from the axis, so no sampled axis division is introduced. Off
support the supplied base third-mixed jet is returned exactly and higher eta
repair data are not requested.

The unrepaired third-mixed ``phi_n`` jet, normalization ``C``, and moment/patch
eta-jets remain upstream inputs. This is Section-5 ``formal-structure`` solver
infrastructure, not a paper-exact recursive coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_phi_jets import _bump_over_radius_x_second_jet
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_preceding_diffusion_parameter_jet import ProfileThirdMixedJet


BasePhiThirdMixedJetProvider = Callable[[float, float], ProfileThirdMixedJet]


def _provider_jet(
    provider: BasePhiThirdMixedJetProvider,
    X: float,
    eta: float,
) -> ProfileThirdMixedJet:
    value = provider(X, eta)
    if not isinstance(value, ProfileThirdMixedJet):
        raise TypeError("base_phi_third_mixed_jet must return ProfileThirdMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedPhiThirdMixedJetAdapter:
    """Lift the compact repaired ``phi_n`` profile to a third-mixed scalar jet."""

    profile: Lemma52RepairedProfileAdapter
    C: float
    base_phi_third_mixed_jet: BasePhiThirdMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired phi=C E/R third-mixed bridge from caller-supplied "
        "unrepaired phi data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        C = float(self.C)
        if not math.isfinite(C) or C <= 0.0:
            raise ValueError("C must be finite and positive")
        if not callable(self.base_phi_third_mixed_jet):
            raise TypeError("base_phi_third_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")
        object.__setattr__(self, "C", C)

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    def phi_third_mixed_jet(self, X: float, eta: float) -> ProfileThirdMixedJet:
        """Return the repaired PositiveAxis ``phi_n`` third-mixed jet."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_phi_third_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.e_bumps)
        if not self._has_active_bump(X, bumps):
            return base

        coeff = self.profile.repair_coefficients(eta, max_order=3)
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
                    "compact phi repair third-mixed correction is outside binary64 range"
                )
            return value

        return ProfileThirdMixedJet(
            value=base.value + correction(0, 0),
            radial=base.radial + correction(0, 1),
            radial2=base.radial2 + correction(0, 2),
            parameter=base.parameter + correction(1, 0),
            radial_parameter=base.radial_parameter + correction(1, 1),
            parameter2=base.parameter2 + correction(2, 0),
            radial2_parameter=base.radial2_parameter + correction(1, 2),
            radial_parameter2=base.radial_parameter2 + correction(2, 1),
            parameter3=base.parameter3 + correction(3, 0),
        )
