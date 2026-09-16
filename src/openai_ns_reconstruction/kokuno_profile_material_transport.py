"""Exact clean-room replay of one public NS-profile transport reduction.

Only short mathematical identities and factual provenance from the public Kokuno
reader are used.  The upstream repository currently advertises no reuse
license, so no checker/proof implementation or substantial authored text is
copied here.

This module verifies a finite algebraic seam at supplied exact jets.  It does
not prove that the imported leading profile exists or satisfies its theorem
hypotheses, and it does not establish a paper-exact velocity or full
reconstruction.
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
SOURCE_PROFILE_PAGES = "24-45,144-157"
SOURCE_CHECKER_PATH = "proof_sources/profiles/exact_checks.py"
SOURCE_CHECKER_SHA256 = "a7ee77024a92d06c1b69055a00b10a7d5fd600ca424fe084e2f4399310cf56e7"
SOURCE_RESULT_SHA256 = "5ef2b5f7cee750088143acda8d4a5ee4f0b12f36f0344ab991f5c9333a842e24"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
SOURCE_DOI = "10.5281/zenodo.22678406"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
IMPORTED_PROFILE_EXISTENCE_PROVED_HERE = False


def _q(value: Exact, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


@dataclass(frozen=True)
class TransportPoint:
    h: Fraction
    eta: Fraction
    X: Fraction

    def __init__(self, *, h: Exact, eta: Exact, X: Exact) -> None:
        hq = _q(h, name="h")
        etaq = _q(eta, name="eta")
        xq = _q(X, name="X")
        if not Fraction(0) < hq < Fraction(1, 2):
            raise ValueError("transport identity requires 0 < h < 1/2")
        if abs(etaq) >= 1:
            raise ValueError("transport identity requires |eta| < 1")
        if xq <= 0:
            raise ValueError("W reduction contains 1/X and requires X > 0")
        object.__setattr__(self, "h", hq)
        object.__setattr__(self, "eta", etaq)
        object.__setattr__(self, "X", xq)

    @property
    def D(self) -> Fraction:
        return Fraction(1, 2) - self.h

    @property
    def d(self) -> Fraction:
        return 1 - self.eta * self.eta

    @property
    def L(self) -> Fraction:
        return 1 - 2 * self.h * self.eta * self.eta


@dataclass(frozen=True)
class ScalarJet:
    value: Fraction
    d_eta: Fraction
    d_X: Fraction

    def __init__(self, *, value: Exact, d_eta: Exact, d_X: Exact) -> None:
        object.__setattr__(self, "value", _q(value, name="value"))
        object.__setattr__(self, "d_eta", _q(d_eta, name="d_eta"))
        object.__setattr__(self, "d_X", _q(d_X, name="d_X"))


def radial_primitive_value(
    point: TransportPoint, *, U: Exact, M: Exact, M_eta: Exact
) -> Fraction:
    """Return the regular-profile V0 algebraic value used in the reduction."""

    uq = _q(U, name="U")
    mq = _q(M, name="M")
    meq = _q(M_eta, name="M_eta")
    return (
        2 * point.eta * point.X * uq
        - 2 * point.D * point.eta * mq
        - point.d * meq
    ) / point.L


def transport_W(point: TransportPoint, *, M: Exact, M_eta: Exact) -> Fraction:
    mq = _q(M, name="M")
    meq = _q(M_eta, name="M_eta")
    return 1 - (2 * point.D * point.eta * mq + point.d * meq) / point.X


def transport_Hc(point: TransportPoint, *, U: Exact) -> Fraction:
    uq = _q(U, name="U")
    return point.D * point.eta + point.d * uq


def composed_material_transport_core(
    point: TransportPoint,
    *,
    b: Exact,
    f: ScalarJet,
    U: Exact,
    V0: Exact,
) -> Fraction:
    """Compose time, radial and axial transport before the W/Hc reduction.

    The returned exact rational is the factor multiplying ``q**(b-1)``.
    The three contributions are kept algebraically separate before summation:
    the fixed-coordinate time chain rule, ``u_r partial_r`` with ``u_r=V0/r``,
    and ``u_z partial_z`` with the leading axial profile factor ``U``.
    """

    bq = _q(b, name="b")
    uq = _q(U, name="U")
    vq = _q(V0, name="V0")

    time_numerator = (
        -bq * f.value
        + point.D * point.eta * f.d_eta
        + point.X * f.d_X
    )
    radial_numerator = point.L * vq * f.d_X
    axial_numerator = uq * (
        2 * bq * point.eta * f.value
        + point.d * f.d_eta
        - 2 * point.eta * point.X * f.d_X
    )
    return (time_numerator + radial_numerator + axial_numerator) / point.L


def reduced_material_transport_core(
    point: TransportPoint,
    *,
    b: Exact,
    f: ScalarJet,
    U: Exact,
    M: Exact,
    M_eta: Exact,
) -> Fraction:
    """Return the public W/Hc form, again without the common q power."""

    bq = _q(b, name="b")
    uq = _q(U, name="U")
    W = transport_W(point, M=M, M_eta=M_eta)
    Hc = transport_Hc(point, U=uq)
    numerator = (
        W * point.X * f.d_X
        + Hc * f.d_eta
        - bq * (1 - 2 * point.eta * uq) * f.value
    )
    return numerator / point.L


def material_transport_reduction_residual(
    point: TransportPoint,
    *,
    b: Exact,
    f: ScalarJet,
    U: Exact,
    M: Exact,
    M_eta: Exact,
    V0: Exact | None = None,
) -> Fraction:
    """Exact residual between composed advection and the reduced W/Hc form."""

    vq = (
        radial_primitive_value(point, U=U, M=M, M_eta=M_eta)
        if V0 is None
        else _q(V0, name="V0")
    )
    return composed_material_transport_core(point, b=b, f=f, U=U, V0=vq) - (
        reduced_material_transport_core(
            point, b=b, f=f, U=U, M=M, M_eta=M_eta
        )
    )


def coefficient_residuals(
    point: TransportPoint,
    *,
    b: Exact,
    U: Exact,
    M: Exact,
    M_eta: Exact,
    V0: Exact | None = None,
) -> tuple[Fraction, Fraction, Fraction]:
    """Return residuals for the f_X, f_eta and f coefficients separately."""

    bq = _q(b, name="b")
    uq = _q(U, name="U")
    vq = (
        radial_primitive_value(point, U=uq, M=M, M_eta=M_eta)
        if V0 is None
        else _q(V0, name="V0")
    )
    W = transport_W(point, M=M, M_eta=M_eta)
    Hc = transport_Hc(point, U=uq)

    x_raw = point.X + point.L * vq - 2 * point.eta * point.X * uq
    eta_raw = point.D * point.eta + point.d * uq
    value_raw = -bq + 2 * bq * point.eta * uq
    return (
        x_raw - W * point.X,
        eta_raw - Hc,
        value_raw + bq * (1 - 2 * point.eta * uq),
    )
