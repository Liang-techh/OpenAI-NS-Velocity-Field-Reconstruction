"""Clean-room exact replay of the public NS-profile radial primitive identity.

The public Kokuno reader displays, with A=1/2+h, D=1/2-h,
d=1-eta^2 and L=1-2h eta^2,

    V0 = (2 eta X U - 2 D eta M - d M_eta) / L,

where M(X,eta)=int_0^X U(x,eta) dx, together with

    d_X V0 = (2 A eta U - d U_eta + 2 eta X U_X) / L.

This module independently reimplements only that short algebra with exact
``Fraction`` inputs. It neither constructs nor proves existence of the imported
leading profile. A small bridge helper compares the exact result with the
target ``LeadingProfile`` API as a floating regression diagnostic only.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, TypeAlias

Exact: TypeAlias = int | Fraction

SOURCE_REPOSITORY = "https://github.com/KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_PUBLIC_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_PROFILE_ARCHIVE_PATH = "proof_sources/profiles/profiles_body.tex"
SOURCE_PROFILE_SHA256 = "63327f4a6d339d39de096230af1810cb83e7c427ba0b0cc78c576b2750437a95"
SOURCE_CHECKER_PATH = "proof_sources/profiles/exact_checks.py"
SOURCE_CHECKER_SHA256 = "a7ee77024a92d06c1b69055a00b10a7d5fd600ca424fe084e2f4399310cf56e7"
SOURCE_RESULT_PATH = "proof_sources/profiles/exact_checks.json"
SOURCE_RESULT_SHA256 = "5ef2b5f7cee750088143acda8d4a5ee4f0b12f36f0344ab991f5c9333a842e24"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
IMPORTED_PROFILE_EXISTENCE_PROVED_HERE = False


def _exact(value: Exact, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


@dataclass(frozen=True)
class RadialPrimitiveJet:
    """Exact point/primitive data sufficient for the displayed V0 identity."""

    X: Fraction
    eta: Fraction
    h: Fraction
    U: Fraction
    U_X: Fraction
    U_eta: Fraction
    M: Fraction
    M_eta: Fraction
    M_X: Fraction
    M_eta_X: Fraction

    def __init__(
        self,
        *,
        X: Exact,
        eta: Exact,
        h: Exact,
        U: Exact,
        U_X: Exact,
        U_eta: Exact,
        M: Exact,
        M_eta: Exact,
        M_X: Exact,
        M_eta_X: Exact,
    ) -> None:
        vals = {
            "X": X,
            "eta": eta,
            "h": h,
            "U": U,
            "U_X": U_X,
            "U_eta": U_eta,
            "M": M,
            "M_eta": M_eta,
            "M_X": M_X,
            "M_eta_X": M_eta_X,
        }
        for name, value in vals.items():
            object.__setattr__(self, name, _exact(value, name=name))
        if self.X < 0:
            raise ValueError("X must be nonnegative")
        if abs(self.eta) > 1:
            raise ValueError("eta must satisfy |eta|<=1")
        if not 0 < self.h < Fraction(1, 2):
            raise ValueError("h must satisfy 0<h<1/2")

    @property
    def A(self) -> Fraction:
        return Fraction(1, 2) + self.h

    @property
    def D(self) -> Fraction:
        return Fraction(1, 2) - self.h

    @property
    def d(self) -> Fraction:
        return 1 - self.eta * self.eta

    @property
    def L(self) -> Fraction:
        return 1 - 2 * self.h * self.eta * self.eta

    @property
    def primitive_compatibility_residuals(self) -> tuple[Fraction, Fraction]:
        """Return M_X-U and (M_eta)_X-U_eta exactly."""
        return self.M_X - self.U, self.M_eta_X - self.U_eta

    def V0(self) -> Fraction:
        numerator = (
            2 * self.eta * self.X * self.U
            - 2 * self.D * self.eta * self.M
            - self.d * self.M_eta
        )
        return numerator / self.L

    def dV0_dX_from_primitive(self) -> Fraction:
        numerator = (
            2 * self.eta * (self.U + self.X * self.U_X)
            - 2 * self.D * self.eta * self.M_X
            - self.d * self.M_eta_X
        )
        return numerator / self.L

    def divergence_rhs(self) -> Fraction:
        numerator = (
            2 * self.A * self.eta * self.U
            - self.d * self.U_eta
            + 2 * self.eta * self.X * self.U_X
        )
        return numerator / self.L

    def divergence_residual(self) -> Fraction:
        return self.dV0_dX_from_primitive() - self.divergence_rhs()


@dataclass(frozen=True)
class LeadingProfileBridgeResult:
    expected_V0: Fraction
    expected_flux_factor: Fraction
    actual_V0: float
    actual_flux_factor: float

    @property
    def V0_error(self) -> float:
        return self.actual_V0 - float(self.expected_V0)

    @property
    def flux_error(self) -> float:
        return self.actual_flux_factor - float(self.expected_flux_factor)


def compare_target_leading_profile(profile: Any, jet: RadialPrimitiveJet) -> LeadingProfileBridgeResult:
    """Compare exact V0 data with the target LeadingProfile API.

    This is intentionally only a floating regression bridge. Exact validation is
    the Fraction identity above; this helper does not upgrade ``profile`` to a
    paper-exact witness.
    """
    if not isinstance(jet, RadialPrimitiveJet):
        raise TypeError("jet must be a RadialPrimitiveJet")
    if jet.X <= 0:
        raise ValueError("target flux-factor bridge requires X>0")
    for attr in ("V0", "radial_flux_factor"):
        if not callable(getattr(profile, attr, None)):
            raise TypeError("profile must implement target LeadingProfile V0/radial_flux_factor")
    provenance = getattr(profile, "provenance", None)
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("profile bridge requires explicit provenance")

    actual_v0 = float(
        profile.V0(
            float(jet.X),
            float(jet.eta),
            float(jet.h),
            d=float(jet.d),
            L=float(jet.L),
        )
    )
    actual_flux = float(
        profile.radial_flux_factor(
            float(jet.X),
            float(jet.eta),
            float(jet.h),
            d=float(jet.d),
            L=float(jet.L),
        )
    )
    expected_v0 = jet.V0()
    return LeadingProfileBridgeResult(
        expected_V0=expected_v0,
        expected_flux_factor=expected_v0 / jet.X,
        actual_V0=actual_v0,
        actual_flux_factor=actual_flux,
    )
