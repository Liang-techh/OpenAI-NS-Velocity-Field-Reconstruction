"""Analytic third eta derivative of the regular Eq. (5.6) ``Omega_k / X`` row.

The landed second-eta implementation remains authoritative for value, first,
and second eta derivatives.  This module differentiates the same five Eq. (5.6)
rows once more using the hierarchy-owned fifth-mixed regular-flux jet introduced
for Stage 2.  Quotients by ``L=1-2h eta^2`` are differentiated by the exact
identity ``L q = N``; production uses no finite differences or fitted data.

This is fail-closed ``formal-structure`` infrastructure.  It does not by itself
materialize the paper's recursive coefficients or establish all-order bounds.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Sequence
import math

from .background_lower_history_source import ProfileSecondJet
from .background_omega_second_parameter_jet import (
    OmegaSecondParameterJet,
    omega_over_x_second_parameter_jet_eq_5_6,
)
from .background_regular_flux_fifth_mixed_jets import RegularFluxFifthMixedJet
from .background_regular_flux_second_jets import AxialThirdMixedJet
from .coordinates import validate_h


@dataclass(frozen=True)
class OmegaThirdParameterJet:
    """Value and first three analytic eta derivatives of ``Omega_k/X``."""

    value: float
    parameter: float
    parameter2: float
    parameter3: float

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2", "parameter3"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


def _order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError("order must be a nonnegative integer")
    return int(value)


def _history(values: Sequence[object], order: int, expected_type: type, name: str) -> tuple[object, ...]:
    if len(values) < order + 1:
        raise ValueError(f"{name} must contain orders 0 through {order}")
    out = tuple(values[j] for j in range(order + 1))
    if any(not isinstance(value, expected_type) for value in out):
        raise TypeError(f"{name} entries must be {expected_type.__name__} instances")
    return out


def _beta_rows(jet: RegularFluxFifthMixedJet) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    eta = (
        jet.value, jet.parameter, jet.parameter2, jet.parameter3,
        jet.parameter4, jet.parameter5,
    )
    radial = (
        jet.radial, jet.radial_parameter, jet.radial_parameter2,
        jet.radial_parameter3, jet.radial_parameter4,
    )
    radial2 = (
        jet.radial2, jet.radial2_parameter, jet.radial2_parameter2,
        jet.radial2_parameter3,
    )
    return eta, radial, radial2


def _quotient_derivatives(
    numerators: Sequence[float], *, h: float, eta: float
) -> tuple[float, ...]:
    """Differentiate ``q=N/L`` exactly from ``L q=N``.

    ``L`` is quadratic in eta, so the recurrence involves only ``L'`` and
    ``L''`` and is exact at every requested finite derivative order.
    """
    h = validate_h(h)
    eta = float(eta)
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("similarity factor L must be positive")
    ell1 = -4.0 * h * eta
    ell2 = -4.0 * h
    out: list[float] = []
    for m, raw in enumerate(numerators):
        rhs = float(raw)
        if m >= 1:
            rhs -= m * ell1 * out[m - 1]
        if m >= 2:
            rhs -= (m * (m - 1) / 2.0) * ell2 * out[m - 2]
        value = rhs / ell
        if not math.isfinite(value):
            raise OverflowError("analytic L-quotient derivative is outside binary64 range")
        out.append(value)
    return tuple(out)


def _z_numerator_derivatives(
    values: Sequence[float],
    radial_values: Sequence[float],
    X: float,
    eta: float,
    *,
    b: float,
    max_order: int,
) -> tuple[float, ...]:
    """Eta derivatives of ``2 eta(b-1)f + d f_eta - 2 eta X f_X``."""
    d = 1.0 - eta * eta
    out: list[float] = []
    for q in range(max_order + 1):
        value = d * values[q + 1] + 2.0 * eta * (b - q - 1.0) * values[q]
        if q:
            value += q * (2.0 * b - q - 1.0) * values[q - 1]
        value -= 2.0 * X * (
            eta * radial_values[q]
            + (q * radial_values[q - 1] if q else 0.0)
        )
        out.append(value)
    return tuple(out)


def _regular_T_parameter3(
    jet: RegularFluxFifthMixedJet, X: float, eta: float, *, b: float, h: float
) -> float:
    D = 0.5 - validate_h(h)
    values, radial, _radial2 = _beta_rows(jet)
    numerators = tuple(
        (1.0 - b + q * D) * values[q]
        + D * eta * values[q + 1]
        + X * radial[q]
        for q in range(4)
    )
    return _quotient_derivatives(numerators, h=h, eta=eta)[3]


def _regular_Z_derivatives(
    jet: RegularFluxFifthMixedJet, X: float, eta: float, *, b: float, h: float
) -> tuple[float, float, float, float]:
    values, radial, _radial2 = _beta_rows(jet)
    numerators = _z_numerator_derivatives(
        values, radial, X, eta, b=b, max_order=3
    )
    result = _quotient_derivatives(numerators, h=h, eta=eta)
    return result[0], result[1], result[2], result[3]


def _regular_Z2_parameter3(
    jet: RegularFluxFifthMixedJet, X: float, eta: float, *, b: float, h: float
) -> float:
    """Return ``partial_eta^3[Z_(b-D) Z_b(X beta) / X]`` analytically."""
    h = validate_h(h)
    D = 0.5 - h
    values, radial, radial2 = _beta_rows(jet)

    # The inner shifted derivative is needed through eta^4; its X derivative is
    # needed through eta^3.  The fifth-mixed beta jet is exactly sufficient.
    inner_num = _z_numerator_derivatives(
        values, radial, X, eta, b=b, max_order=4
    )
    inner_x_num = _z_numerator_derivatives(
        radial, radial2, X, eta, b=b - 1.0, max_order=3
    )
    inner = _quotient_derivatives(inner_num, h=h, eta=eta)
    inner_x = _quotient_derivatives(inner_x_num, h=h, eta=eta)

    outer_num = _z_numerator_derivatives(
        inner, inner_x, X, eta, b=b - D, max_order=3
    )
    return _quotient_derivatives(outer_num, h=h, eta=eta)[3]


def _axial_second(jet: AxialThirdMixedJet) -> ProfileSecondJet:
    return ProfileSecondJet(
        value=jet.value,
        radial=jet.radial,
        radial2=jet.radial2,
        parameter=jet.parameter,
        radial_parameter=jet.radial_parameter,
        parameter2=jet.parameter2,
    )


def omega_over_x_third_parameter_jet_eq_5_6(
    order: int,
    X: float,
    eta: float,
    radial_flux_jets: Sequence[RegularFluxFifthMixedJet],
    axial_jets: Sequence[AxialThirdMixedJet],
    *,
    h: float,
) -> OmegaThirdParameterJet:
    """Return value through ``partial_eta^3(Omega_order/X)`` from Eq. (5.6)."""
    order = _order(order)
    beta = _history(radial_flux_jets, order, RegularFluxFifthMixedJet, "radial_flux_jets")
    axial = _history(axial_jets, order, AxialThirdMixedJet, "axial_jets")
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    h = validate_h(h)

    lower: OmegaSecondParameterJet = omega_over_x_second_parameter_jet_eq_5_6(
        order,
        X,
        eta,
        [jet.fourth() for jet in beta],
        [_axial_second(jet) for jet in axial],
        h=h,
    )
    lam = lambda n: 2.0 * n * h
    parameter3 = _regular_T_parameter3(beta[order], X, eta, b=lam(order), h=h)

    for i in range(order + 1):
        j = order - i
        left, _left_x, _left_xx = _beta_rows(beta[i])
        right_beta, right_x, _right_xx = _beta_rows(beta[j])
        right = tuple(0.5 * right_beta[q] + X * right_x[q] for q in range(4))
        parameter3 += (
            left[3] * right[0]
            + 3.0 * left[2] * right[1]
            + 3.0 * left[1] * right[2]
            + left[0] * right[3]
        )

        z = _regular_Z_derivatives(beta[j], X, eta, b=lam(j), h=h)
        a = (axial[i].value, axial[i].parameter, axial[i].parameter2, axial[i].parameter3)
        parameter3 += a[3] * z[0] + 3.0 * a[2] * z[1] + 3.0 * a[1] * z[2] + a[0] * z[3]

    _v, x, xx = _beta_rows(beta[order])
    parameter3 += -2.0 * (2.0 * x[3] + X * xx[3])
    if order >= 1:
        parameter3 -= _regular_Z2_parameter3(
            beta[order - 1], X, eta, b=lam(order - 1), h=h
        )

    if not math.isfinite(parameter3):
        raise OverflowError("Eq. (5.6) omega third parameter jet is outside binary64 range")
    return OmegaThirdParameterJet(
        value=lower.value,
        parameter=lower.parameter,
        parameter2=lower.parameter2,
        parameter3=parameter3,
    )
