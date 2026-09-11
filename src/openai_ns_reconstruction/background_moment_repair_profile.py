"""Function-level bridge from Lemma 5.2 repair jets to repaired profiles.

The pointwise repair in :mod:`background_moment_repair` and its eta-jet lift in
:mod:`background_moment_repair_jets` already solve Eqs. (5.14)-(5.16) once the
five unrepaired moments and the patch factor ``p(eta)=e_* f(eta)`` are supplied.
This module turns those coefficient jets into actual repaired ``U_n``/``E_n``
functions, including the eta derivative of ``U_n`` required by Eq. (5.2) and
Eq. (5.15).

The adapter deliberately keeps the upstream moment jets and patch-factor jets as
explicit caller inputs.  Until the materialized Section-4/5 hierarchy supplies
those functions with uniform analytic/nonvanishing certificates, the resulting
profile is formal-structure infrastructure and is never marked paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

import numpy as np

from .background_moment_repair import Lemma52MomentRepair
from .background_moment_repair_jets import (
    Lemma52RepairJetCoefficients,
    solve_lemma52_eta_jets,
)
from .profiles import LeadingProfile, ScalarFn

MomentJetProvider = Callable[[float], object]
PatchJetProvider = Callable[[float], object]


def _eta(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or abs(value) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return value


def _order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("max_order must be a nonnegative integer")
    return value


@dataclass(frozen=True)
class Lemma52RepairedProfileAdapter:
    """Apply the compact five-bump repair as eta-dependent profile functions.

    ``base_moment_derivatives(eta)`` must return an array with shape
    ``(J+1, 5)`` containing ordinary eta derivatives of the five unrepaired
    moments. ``patch_factor_derivatives(eta)`` must return ``(J+1,)`` ordinary
    derivatives of ``p(eta)=e_* f(eta)``.  Only the rows requested by
    :meth:`repair_coefficients` are consumed.

    The radial bump functions are independent of eta, so after differentiating
    Eq. (5.14),

        d_eta U_n = d_eta U_tilde_n + sum_j alpha'_j(eta) b^U_j(R).

    This is the derivative needed by the landed Eq. (5.2)/(5.15) reconstruction.
    """

    repair: Lemma52MomentRepair
    base_profile: LeadingProfile
    base_moment_derivatives: MomentJetProvider
    patch_factor_derivatives: PatchJetProvider
    provenance: str = (
        "Lemma 5.2 functional repair adapter from caller-supplied moment/patch "
        "jets; formal-structure only, not a paper-exact coefficient hierarchy."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.repair, Lemma52MomentRepair):
            raise TypeError("repair must be a Lemma52MomentRepair")
        if not isinstance(self.base_profile, LeadingProfile):
            raise TypeError("base_profile must be a LeadingProfile")
        if not callable(self.base_moment_derivatives):
            raise TypeError("base_moment_derivatives must be callable")
        if not callable(self.patch_factor_derivatives):
            raise TypeError("patch_factor_derivatives must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

    def repair_coefficients(
        self,
        eta: float,
        *,
        max_order: int = 1,
    ) -> Lemma52RepairJetCoefficients:
        """Return repair coefficient derivatives through ``max_order`` at eta."""

        eta = _eta(eta)
        max_order = _order(max_order)
        moments = np.asarray(self.base_moment_derivatives(eta), dtype=float)
        patch = np.asarray(self.patch_factor_derivatives(eta), dtype=float)
        need = max_order + 1
        if moments.ndim != 2 or moments.shape[1:] != (5,) or moments.shape[0] < need:
            raise ValueError(
                "base_moment_derivatives must provide at least max_order+1 rows of 5 moments"
            )
        if patch.ndim != 1 or patch.shape[0] < need:
            raise ValueError(
                "patch_factor_derivatives must provide at least max_order+1 entries"
            )
        if not np.all(np.isfinite(moments[:need])) or not np.all(np.isfinite(patch[:need])):
            raise ValueError("requested moment and patch-factor jets must be finite")
        return solve_lemma52_eta_jets(
            self.repair,
            moments[:need],
            patch[:need],
        )

    def U(self, X: float, eta: float) -> float:
        """Evaluate the repaired axial coefficient ``U_n`` from Eq. (5.14)."""

        X, eta = self.base_profile._point(X, eta)
        coeff = self.repair_coefficients(eta, max_order=0)
        radius = math.sqrt(2.0 * X)
        correction = math.fsum(
            float(a) * bump(radius)
            for a, bump in zip(coeff.alpha_derivatives[0], self.repair.u_bumps)
        )
        value = float(self.base_profile.U(X, eta)) + correction
        if not math.isfinite(value):
            raise OverflowError("repaired U_n is outside floating-point range")
        return value

    def dU_deta(self, X: float, eta: float) -> float:
        """Evaluate ``d_eta U_n`` using the differentiated compact repair."""

        X, eta = self.base_profile._point(X, eta)
        coeff = self.repair_coefficients(eta, max_order=1)
        radius = math.sqrt(2.0 * X)
        correction = math.fsum(
            float(a) * bump(radius)
            for a, bump in zip(coeff.alpha_derivatives[1], self.repair.u_bumps)
        )
        value = float(self.base_profile.dU_deta(X, eta)) + correction
        if not math.isfinite(value):
            raise OverflowError("repaired d_eta U_n is outside floating-point range")
        return value

    def E(self, X: float, eta: float) -> float:
        """Evaluate the repaired swirl coefficient ``E_n`` from Eq. (5.14)."""

        X, eta = self.base_profile._point(X, eta)
        coeff = self.repair_coefficients(eta, max_order=0)
        radius = math.sqrt(2.0 * X)
        correction = math.fsum(
            float(beta) * bump(radius)
            for beta, bump in zip(coeff.beta_derivatives[0], self.repair.e_bumps)
        )
        value = float(self.base_profile.E(X, eta)) + correction
        if not math.isfinite(value):
            raise OverflowError("repaired E_n is outside floating-point range")
        return value

    def F(self, X: float, eta: float) -> float:
        """Return the axis-regular swirl factor ``E_n/sqrt(2X)``.

        The compact repair bumps have strictly positive radial support, so the
        axis value is unchanged and may be inherited from the base profile's
        certified/supplied smooth ``F``.  Away from the axis the corrected value
        is evaluated directly from the repaired ``E_n``.
        """

        X, eta = self.base_profile._point(X, eta)
        if X == 0.0:
            return self.base_profile.smooth_swirl_factor(0.0, eta)
        value = self.E(X, eta) / math.sqrt(2.0 * X)
        if not math.isfinite(value):
            raise OverflowError("repaired smooth swirl factor is nonfinite")
        return value

    def as_leading_profile(
        self,
        *,
        Pi: ScalarFn | None = None,
        name: str | None = None,
    ) -> LeadingProfile:
        """Expose the repaired coefficient through the repository profile API.

        ``paper_exact`` is intentionally hard-coded to ``False``.  ``Pi`` is an
        optional independently reconstructed/supplied pressure profile; it is
        never inherited silently from the unrepaired base profile.
        """

        if Pi is not None and not callable(Pi):
            raise TypeError("Pi must be callable or None")
        return LeadingProfile(
            E=self.E,
            U=self.U,
            dU_deta=self.dU_deta,
            Pi=Pi,
            name=name or f"lemma52-repaired:{self.base_profile.name}",
            paper_exact=False,
            F=self.F if self.base_profile.F is not None else None,
            provenance=self.provenance,
        )
