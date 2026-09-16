"""Independent exact replay of one public Kokuno base-heat identity.

The Kokuno workbench does not currently expose a clear repository reuse license.
This module therefore does not copy ``verify_identities.py``.  It independently
reimplements only the polynomial reduction printed in the public workbench for
the base-heat factor ODE, using exact ``Fraction`` coefficient arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Mapping, TypeAlias

Monomial: TypeAlias = tuple[int, int, int]  # powers of (h, Z, v)
Polynomial: TypeAlias = Mapping[Monomial, Fraction]

SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_WORKBENCH_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_WORKBENCH_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
BASE_HEAT_DERIVATION_SHA256 = "98d101542a7d8d5c1475764ec99178a663fe363e416d572eeafa048e724a193d"
BASE_HEAT_CHECKER_SHA256 = "da6fe495af223772824b2f658c83c4aec08d7aa41882818589eaf14c6a3425dc"
BASE_HEAT_RECORDED_RESULT_SHA256 = "5b74d2d2598841ce2512790dd70efe35640b69a4ea5a8516a4dbbd482f0bd843"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_WORKBENCH_URL = (
    "https://github.com/KokunoYumeto/yang-mills-interacting-workbench/"
    f"blob/{SOURCE_WORKBENCH_COMMIT}/navier-stokes/navier_stokes_workbench.tex"
)

_ZERO_MONOMIAL: Monomial = (0, 0, 0)


def _normalize(poly: Mapping[Monomial, Fraction]) -> dict[Monomial, Fraction]:
    return {key: Fraction(value) for key, value in poly.items() if value != 0}


def _constant(value: int | Fraction) -> dict[Monomial, Fraction]:
    coefficient = Fraction(value)
    return {} if coefficient == 0 else {_ZERO_MONOMIAL: coefficient}


def _variable(index: int) -> dict[Monomial, Fraction]:
    if index not in (0, 1, 2):
        raise ValueError("polynomial variable index must be 0, 1, or 2")
    exponent = [0, 0, 0]
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def _add(*polys: Mapping[Monomial, Fraction]) -> dict[Monomial, Fraction]:
    out: dict[Monomial, Fraction] = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + Fraction(coefficient)
            if out[monomial] == 0:
                del out[monomial]
    return out


def _scale(poly: Mapping[Monomial, Fraction], scalar: int | Fraction) -> dict[Monomial, Fraction]:
    scalar_fraction = Fraction(scalar)
    return _normalize({key: value * scalar_fraction for key, value in poly.items()})


def _mul(*polys: Mapping[Monomial, Fraction]) -> dict[Monomial, Fraction]:
    out = _constant(1)
    for poly in polys:
        product: dict[Monomial, Fraction] = {}
        for left_power, left_coefficient in out.items():
            for right_power, right_coefficient in poly.items():
                exponent = tuple(left_power[i] + right_power[i] for i in range(3))
                product[exponent] = product.get(exponent, Fraction(0)) + left_coefficient * right_coefficient
        out = _normalize(product)
    return out


def _square(poly: Mapping[Monomial, Fraction]) -> dict[Monomial, Fraction]:
    return _mul(poly, poly)


def same_polynomial(left: Mapping[Monomial, Fraction], right: Mapping[Monomial, Fraction]) -> bool:
    """Return exact coefficient equality after removing zero coefficients."""

    return _normalize(left) == _normalize(right)


def _freeze(poly: Mapping[Monomial, Fraction]) -> Polynomial:
    return MappingProxyType(_normalize(poly))


@dataclass(frozen=True)
class BaseHeatOdeReplay:
    """Exact outputs for the independently reconstructed ODE sub-check."""

    ode_common_factor_polynomial: Polynomial
    reduced_polynomial: Polynomial
    total_derivative_factor_polynomial: Polynomial
    radial_heat_coefficient_lhs: Polynomial
    radial_heat_coefficient_rhs: Polynomial
    tolerance: Fraction

    @property
    def ode_reduction_exact(self) -> bool:
        return same_polynomial(self.ode_common_factor_polynomial, self.reduced_polynomial)

    @property
    def total_derivative_factor_exact(self) -> bool:
        return same_polynomial(self.total_derivative_factor_polynomial, self.reduced_polynomial)

    @property
    def radial_heat_coefficient_exact(self) -> bool:
        return same_polynomial(self.radial_heat_coefficient_lhs, self.radial_heat_coefficient_rhs)

    @property
    def passed(self) -> bool:
        return (
            self.tolerance == 0
            and self.ode_reduction_exact
            and self.total_derivative_factor_exact
            and self.radial_heat_coefficient_exact
        )


def replay_base_heat_ode_polynomial() -> BaseHeatOdeReplay:
    """Replay the public base-heat ODE reduction as a universal polynomial identity.

    The public workbench writes, after factoring the ODE integrand,

        (1+h) Z^2 v^2
        - (1 + 2(1+h)Z) v (1+Zv)
        + (1+h)(1+Zv)^2
        = (1+h) - v - Z v^2.

    It also identifies the right-hand side as the polynomial left after
    differentiating ``exp(-v) v^(1+h) (1+Zv)^(-1-h)`` by ``v`` with the
    same common factor.  This implementation expands both statements over
    Q[h,Z,v] exactly.  No numerical quadrature or tolerance is used.
    """

    one = _constant(1)
    half = Fraction(1, 2)
    h = _variable(0)
    z = _variable(1)
    v = _variable(2)
    one_plus_h = _add(one, h)
    one_plus_zv = _add(one, _mul(z, v))

    ode_polynomial = _add(
        _mul(one_plus_h, _square(z), _square(v)),
        _scale(
            _mul(
                _add(one, _scale(_mul(one_plus_h, z), 2)),
                v,
                one_plus_zv,
            ),
            -1,
        ),
        _mul(one_plus_h, _square(one_plus_zv)),
    )

    reduced = _add(
        one_plus_h,
        _scale(v, -1),
        _scale(_mul(z, _square(v)), -1),
    )

    derivative_factor = _add(
        _scale(_mul(v, one_plus_zv), -1),
        _mul(one_plus_h, one_plus_zv),
        _scale(_mul(one_plus_h, z, v), -1),
    )

    a = _add(_constant(half), h)
    radial_lhs = _add(_scale(_square(a), 2), _constant(-half))
    radial_rhs = _scale(_mul(h, one_plus_h), 2)

    return BaseHeatOdeReplay(
        ode_common_factor_polynomial=_freeze(ode_polynomial),
        reduced_polynomial=_freeze(reduced),
        total_derivative_factor_polynomial=_freeze(derivative_factor),
        radial_heat_coefficient_lhs=_freeze(radial_lhs),
        radial_heat_coefficient_rhs=_freeze(radial_rhs),
        tolerance=Fraction(0),
    )
