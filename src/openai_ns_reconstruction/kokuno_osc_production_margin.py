"""Clean-room exact replay of the corrected oscillatory production error budget.

This module independently implements only the short public AS35/AS37 arithmetic
from the corrected Kokuno reader.  It does not copy or execute the source
checker and it does not construct the oscillatory pulse.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


SOURCE_H_MAX = Fraction(1, 100)


def _fraction(name: str, value: Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


def _positive(name: str, value: Fraction) -> Fraction:
    value = _fraction(name, value)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _nonnegative(name: str, value: Fraction) -> Fraction:
    value = _fraction(name, value)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _ell(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("ell must be an integer")
    if value < 1:
        raise ValueError("ell must be positive")
    return value


@dataclass(frozen=True)
class OscillatoryProductionBudget:
    """Exact AS35 error budget and the resulting AS37 lower coefficient."""

    h: Fraction
    ell: int
    dyadic_exponent: int
    frozen_frame_error: Fraction
    borel_actual_shear_error: Fraction
    box_actual_shear_error: Fraction
    stress_factor: Fraction
    actual_shear_error: Fraction
    total_error: Fraction
    half_leading_budget: Fraction
    half_leading_margin: Fraction
    lower_production_coefficient: Fraction

    @property
    def certifies_half_leading_bound(self) -> bool:
        return self.half_leading_margin >= 0


def replay_oscillatory_production_budget(
    *,
    h: Fraction,
    ell: int,
    g_plus: Fraction,
    g_minus: Fraction,
    c_minus: Fraction,
    c_t: Fraction,
    c_b: Fraction,
    l_g: Fraction,
    c_box: Fraction,
    s_max: Fraction,
    gamma_max: Fraction,
    c_x: Fraction,
) -> OscillatoryProductionBudget:
    """Evaluate the corrected production budget exactly on a dyadic lattice.

    The public corrected reader uses ``S_* = ell^2`` and bounds the AS35 error
    by

      g_+ C_T ell^-2
      + (C_B 2^(-2 h ell) + L_g C_box ell^-6)
        (s_max + Gamma_max + C_T ell^-2).

    Exact rational arithmetic can represent the dyadic factor without an
    approximation when ``2*h*ell`` is an integer.  This focused replay rejects
    other points rather than silently approximating them.
    """

    h = _positive("h", h)
    if h >= SOURCE_H_MAX:
        raise ValueError("source hypothesis requires 0 < h < 1/100")
    ell = _ell(ell)
    g_plus = _nonnegative("g_plus", g_plus)
    g_minus = _positive("g_minus", g_minus)
    c_minus = _positive("c_minus", c_minus)
    c_t = _nonnegative("c_t", c_t)
    c_b = _nonnegative("c_b", c_b)
    l_g = _nonnegative("l_g", l_g)
    c_box = _nonnegative("c_box", c_box)
    s_max = _nonnegative("s_max", s_max)
    gamma_max = _nonnegative("gamma_max", gamma_max)
    c_x = _positive("c_x", c_x)

    exponent = 2 * h * ell
    if exponent.denominator != 1:
        raise ValueError("exact replay requires integral 2*h*ell")
    dyadic_exponent = exponent.numerator
    if dyadic_exponent < 0:
        raise ValueError("dyadic exponent must be nonnegative")

    ell2 = Fraction(ell * ell, 1)
    ell6 = Fraction(ell**6, 1)
    frozen_frame_error = g_plus * c_t / ell2
    borel_error = c_b * Fraction(1, 2**dyadic_exponent)
    box_error = l_g * c_box / ell6
    stress_factor = s_max + gamma_max + c_t / ell2
    actual_shear_error = (borel_error + box_error) * stress_factor
    total_error = frozen_frame_error + actual_shear_error

    half_leading_budget = g_minus * c_minus / 2
    half_leading_margin = half_leading_budget - total_error
    lower_production_coefficient = half_leading_budget * c_x * c_x

    return OscillatoryProductionBudget(
        h=h,
        ell=ell,
        dyadic_exponent=dyadic_exponent,
        frozen_frame_error=frozen_frame_error,
        borel_actual_shear_error=borel_error,
        box_actual_shear_error=box_error,
        stress_factor=stress_factor,
        actual_shear_error=actual_shear_error,
        total_error=total_error,
        half_leading_budget=half_leading_budget,
        half_leading_margin=half_leading_margin,
        lower_production_coefficient=lower_production_coefficient,
    )


def require_half_leading_production(**kwargs: object) -> OscillatoryProductionBudget:
    """Fail closed unless the full corrected AS35 budget meets AS37's half bound."""

    result = replay_oscillatory_production_budget(**kwargs)  # type: ignore[arg-type]
    if not result.certifies_half_leading_bound:
        raise ValueError(
            "corrected actual-shear error exceeds the half-leading production budget: "
            f"margin={result.half_leading_margin}"
        )
    return result
