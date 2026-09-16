"""Independent exact checks for the public NS-profile coordinate identities.

This module is a clean-room reimplementation from publicly displayed formulas in
KokunoYumeto/yang-mills-interacting-workbench. It does not copy the source
checker or proof body. The source repository does not currently advertise a
reuse license, so only mathematical identities and provenance metadata are
retained here.

Truth boundary: these helpers check finite algebraic coordinate identities.
They do not prove existence, regularity, decay, normalization, admissibility, or
assembly of the imported leading profiles, and they do not establish a
paper-exact velocity field or full reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import TypeAlias

Exact: TypeAlias = int | Fraction

SOURCE_REPOSITORY = "https://github.com/KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_PUBLIC_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_PROFILE_ARCHIVE_PATH = "proof_sources/profiles/profiles_body.tex"
SOURCE_PROFILE_SHA256 = "63327f4a6d339d39de096230af1810cb83e7c427ba0b0cc78c576b2750437a95"
SOURCE_CHECKER_PATH = "proof_sources/profiles/exact_checks.py"
SOURCE_CHECKER_SHA256 = "a7ee77024a92d06c1b69055a00b10a7d5fd600ca424fe084e2f4399310cf56e7"
SOURCE_CHECK_RESULT_PATH = "proof_sources/profiles/exact_checks.json"
SOURCE_CHECK_RESULT_SHA256 = "5ef2b5f7cee750088143acda8d4a5ee4f0b12f36f0344ab991f5c9333a842e24"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
SOURCE_MANUSCRIPT_SHA256 = "8c8a94ad9ac824c8b605b9827cadf7beaca48bd10b380de3cfc872a2c37afa81"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
IMPORTED_PROFILE_EXISTENCE_PROVED_HERE = False


def _fraction(value: Exact, *, name: str) -> Fraction:
    """Convert an exact scalar while rejecting binary floating-point inputs."""

    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


@dataclass(frozen=True)
class SelfSimilarPoint:
    """Interior point of the public ``(q, eta, X)`` coordinate chart.

    The public reader states the coordinate lemma for ``0 < h < 1/2`` and
    physical interior points ``q > 0``, ``|eta| < 1``. Profile endpoint values
    ``eta = +/-1`` are one-sided parameter endpoints and are deliberately not
    accepted as physical chart points here.
    """

    h: Fraction
    q: Fraction
    eta: Fraction
    X: Fraction

    def __init__(self, *, h: Exact, q: Exact, eta: Exact, X: Exact) -> None:
        h_q = _fraction(h, name="h")
        q_q = _fraction(q, name="q")
        eta_q = _fraction(eta, name="eta")
        x_q = _fraction(X, name="X")
        if not Fraction(0) < h_q < Fraction(1, 2):
            raise ValueError("coordinate identity requires 0 < h < 1/2")
        if q_q <= 0:
            raise ValueError("physical coordinate chart requires q > 0")
        if abs(eta_q) >= 1:
            raise ValueError("physical coordinate chart requires |eta| < 1")
        if x_q < 0:
            raise ValueError("axisymmetric coordinate X must be nonnegative")
        object.__setattr__(self, "h", h_q)
        object.__setattr__(self, "q", q_q)
        object.__setattr__(self, "eta", eta_q)
        object.__setattr__(self, "X", x_q)

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
    def scaled_tau_z_jacobian(self) -> Fraction:
        """Return ``det d(tau,z)/d(q,eta) / q**D`` exactly.

        The public map ``tau=q(1-eta**2)``, ``z=q**D eta`` has determinant
        ``q**D * L``. Returning the rational factor ``L`` avoids introducing
        an inexact fractional power of an arbitrary rational ``q``.
        """

        return self.L

    def inverse_derivative_factors(self) -> "InverseDerivativeFactors":
        """Exact rational factors in the public fixed-coordinate derivatives.

        Powers of ``q`` are factored out exactly as displayed in the source:
        ``q_t = q_t_factor``, ``eta_t = eta_t_factor / q``,
        ``X_t = X_t_factor / q``, ``q_z = q**(1-D)*q_z_factor``,
        ``eta_z = q**(-D)*eta_z_factor``, and
        ``X_z = q**(-D)*X_z_factor``.
        """

        L = self.L
        return InverseDerivativeFactors(
            q_t_factor=-1 / L,
            eta_t_factor=self.D * self.eta / L,
            X_t_factor=self.X / L,
            q_z_factor=2 * self.eta / L,
            eta_z_factor=self.d / L,
            X_z_factor=-2 * self.eta * self.X / L,
        )


@dataclass(frozen=True)
class InverseDerivativeFactors:
    q_t_factor: Fraction
    eta_t_factor: Fraction
    X_t_factor: Fraction
    q_z_factor: Fraction
    eta_z_factor: Fraction
    X_z_factor: Fraction


@dataclass(frozen=True)
class ProfileJet:
    """A value and its ``eta``/``X`` derivatives at one profile point."""

    value: Fraction
    d_eta: Fraction
    d_X: Fraction

    def __init__(self, value: Exact, d_eta: Exact, d_X: Exact) -> None:
        object.__setattr__(self, "value", _fraction(value, name="value"))
        object.__setattr__(self, "d_eta", _fraction(d_eta, name="d_eta"))
        object.__setattr__(self, "d_X", _fraction(d_X, name="d_X"))


def scaled_profile_time_derivative(
    point: SelfSimilarPoint, *, b: Exact, jet: ProfileJet
) -> Fraction:
    """Return the exact factor multiplying ``q**(b-1)`` in ``d_t(q**b f)``."""

    b_q = _fraction(b, name="b")
    numerator = (
        -b_q * jet.value
        + point.D * point.eta * jet.d_eta
        + point.X * jet.d_X
    )
    return numerator / point.L


def scaled_profile_axial_derivative(
    point: SelfSimilarPoint, *, b: Exact, jet: ProfileJet
) -> Fraction:
    """Return the exact factor multiplying ``q**(b-D)`` in ``d_z(q**b f)``."""

    b_q = _fraction(b, name="b")
    numerator = (
        2 * b_q * point.eta * jet.value
        + point.d * jet.d_eta
        - 2 * point.eta * point.X * jet.d_X
    )
    return numerator / point.L


def incompressibility_profile_residual(
    point: SelfSimilarPoint,
    *,
    V0_X: Exact,
    U: ProfileJet,
) -> Fraction:
    """Exact residual of the public scalar incompressibility identity.

    Zero is the displayed relation
    ``(V0)_X = L^-1(2 A eta U - d U_eta + 2 eta X U_X)``.
    No profile existence or regularity statement is inferred from a zero
    residual at finitely supplied jets.
    """

    lhs = _fraction(V0_X, name="V0_X")
    rhs = (
        2 * point.A * point.eta * U.value
        - point.d * U.d_eta
        + 2 * point.eta * point.X * U.d_X
    ) / point.L
    return lhs - rhs


@dataclass(frozen=True)
class AxisymmetricFieldJet:
    """Physical axisymmetric derivatives needed for ``omega_theta / r``."""

    s: Fraction
    d_z_ru_r: Fraction
    d_s_u_z: Fraction

    def __init__(self, *, s: Exact, d_z_ru_r: Exact, d_s_u_z: Exact) -> None:
        s_q = _fraction(s, name="s")
        if s_q <= 0:
            raise ValueError("omega_theta/r quotient requires s=r^2/2 > 0")
        object.__setattr__(self, "s", s_q)
        object.__setattr__(self, "d_z_ru_r", _fraction(d_z_ru_r, name="d_z_ru_r"))
        object.__setattr__(self, "d_s_u_z", _fraction(d_s_u_z, name="d_s_u_z"))


def omega_theta_over_r(jet: AxisymmetricFieldJet) -> Fraction:
    """Return ``omega_theta/r`` from ``s=r^2/2`` axisymmetric field data.

    With ``V0=r*u_r`` and no angular dependence,
    ``omega_theta/r = (partial_z V0)/(2s) - partial_s u_z`` exactly.
    This is a component identity only; it does not assert that supplied jets
    arise from an imported leading profile.
    """

    return jet.d_z_ru_r / (2 * jet.s) - jet.d_s_u_z
