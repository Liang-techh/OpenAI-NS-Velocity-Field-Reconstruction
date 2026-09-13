"""Exact third-eta derivative primitive for the Section 5 Eq. (5.7) forcing.

The repaired hierarchy currently owns the exact forcing only through
``partial_eta^2 f_n``.  The next differentiated Picard step requires
``partial_eta^3 f_n``.  This module isolates the purely algebraic forcing part
of that requirement: given an analytic strict-lower source jet through third
eta order, it returns the exact third eta derivative of the displayed Eq. (5.7)
forcing.

This is deliberately *not* called hierarchy-owned.  The third source jet is an
explicit structural input until the Section-5 lower-history hierarchy itself
owns that derivative.  Production uses no finite differences, sampled
coefficient tables, generic cutoff, or fitted data, and this module makes no
paper-exact or Picard-convergence claim.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .background_positive_axis import PositiveAxisSourceJet, pressure_source
from .coordinates import validate_h


@dataclass(frozen=True)
class PositiveAxisSourceThirdParameterJet:
    """Strict-lower source and its first three analytic eta derivatives."""

    value: PositiveAxisSourceJet
    parameter: PositiveAxisSourceJet
    parameter2: PositiveAxisSourceJet
    parameter3: PositiveAxisSourceJet

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2", "parameter3"):
            if not isinstance(getattr(self, name), PositiveAxisSourceJet):
                raise TypeError(f"{name} must be a PositiveAxisSourceJet")


def _positive_C(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("C must be finite and positive")
    return value


def _xi(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    return value


def _eta(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or abs(value) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return value


def positive_axis_forcing_third_parameter_from_source_jet(
    h: float,
    C: float,
    xi: float,
    eta: float,
    source: PositiveAxisSourceThirdParameterJet,
) -> np.ndarray:
    """Return ``partial_eta^3 f_n`` from an analytic third-order source jet.

    For the only nonlinear-in-eta forcing geometry, write
    ``g(eta)=eta/ell`` with ``ell=1-2 h eta^2``.  The exact derivatives used
    here are

    ``g'   = 1/ell + 4 h eta^2/ell^2``,
    ``g''  = 12 h eta/ell^2 + 32 h^2 eta^3/ell^3``,
    ``g''' = 12 h/ell^2 + 192 h^2 eta^2/ell^3
             + 384 h^3 eta^4/ell^4``.

    Hence the sixth row uses the exact third-order Leibniz rule for
    ``g * pressureSource``.  The first three rows remain identically zero.
    """

    h = validate_h(h)
    C = _positive_C(C)
    xi = _xi(xi)
    eta = _eta(eta)
    if not isinstance(source, PositiveAxisSourceThirdParameterJet):
        raise TypeError("source must be a PositiveAxisSourceThirdParameterJet")

    p0 = pressure_source(C, source.value)
    p1 = pressure_source(C, source.parameter)
    p2 = pressure_source(C, source.parameter2)
    p3 = pressure_source(C, source.parameter3)

    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("ell(h,eta) must be positive")
    X = xi * xi

    g = eta / ell
    g1 = 1.0 / ell + 4.0 * h * eta * eta / (ell * ell)
    g2 = (
        12.0 * h * eta / (ell * ell)
        + 32.0 * h * h * eta * eta * eta / (ell * ell * ell)
    )
    g3 = (
        12.0 * h / (ell * ell)
        + 192.0 * h * h * eta * eta / (ell * ell * ell)
        + 384.0
        * h
        * h
        * h
        * eta
        * eta
        * eta
        * eta
        / (ell * ell * ell * ell)
    )

    out = np.zeros(6, dtype=float)
    out[3] = 2.0 * xi * p3
    out[4] = 2.0 * source.parameter3.angular
    out[5] = (
        2.0 * source.parameter3.axial
        - 4.0
        * X
        * (
            g3 * p0
            + 3.0 * g2 * p1
            + 3.0 * g1 * p2
            + g * p3
        )
    )
    if not np.all(np.isfinite(out)):
        raise OverflowError("positive-axis forcing third parameter is outside binary64 range")
    out.setflags(write=False)
    return out
