"""Independent exact replay of the corrected oscillatory actual-shear split.

This module does not copy the Kokuno checker.  It implements only the short
public mathematical identity recorded in the corrected workbench, using exact
rational coordinates in the (N, K) shear plane.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Rational


ExactInput = int | Fraction


def _exact_fraction(value: ExactInput, *, name: str) -> Fraction:
    """Return an exact rational theorem input, rejecting approximate scalars."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    # Deliberately do not accept float/Decimal or generic approximate Rational
    # implementations.  The checker is an exact algebra replay.
    if isinstance(value, Rational):
        return Fraction(value.numerator, value.denominator)
    raise TypeError(f"{name} must be an exact integer/Fraction")


@dataclass(frozen=True)
class ActualShearReplayResult:
    """Exact terms in the corrected actual-shear production decomposition."""

    frozen_shear_magnitude: Fraction
    delta_normal: Fraction
    delta_tangent: Fraction
    flux_normal: Fraction
    flux_tangent: Fraction
    lhs_actual_production: Fraction
    leading_frozen_term: Fraction
    retained_remainder: Fraction
    rhs_corrected_split: Fraction
    identity_residual: Fraction
    leading_only_defect: Fraction

    @property
    def exact_identity_verified(self) -> bool:
        return self.identity_residual == 0

    @property
    def leading_only_is_valid(self) -> bool:
        return self.leading_only_defect == 0


def replay_actual_shear_identity(
    *,
    frozen_shear_magnitude: ExactInput,
    delta_normal: ExactInput,
    delta_tangent: ExactInput,
    flux_normal: ExactInput,
    flux_tangent: ExactInput,
) -> ActualShearReplayResult:
    """Replay ``-g·T = -|g0| T_N - (g-g0)·T`` exactly.

    Coordinates use the public corrected-edition decomposition

    ``g0 = |g0| N``, ``g-g0 = delta_normal N + delta_tangent K``, and
    ``T = flux_normal N + flux_tangent K``.

    Only the Euclidean dot-product algebra in this displayed orthonormal basis
    is checked here.  The routine does not construct the pulse, prove the
    positive-production inequalities, or validate the complete oscillatory
    source component.
    """

    magnitude = _exact_fraction(frozen_shear_magnitude, name="frozen_shear_magnitude")
    if magnitude < 0:
        raise ValueError("frozen_shear_magnitude represents |g0| and must be nonnegative")

    d_n = _exact_fraction(delta_normal, name="delta_normal")
    d_k = _exact_fraction(delta_tangent, name="delta_tangent")
    t_n = _exact_fraction(flux_normal, name="flux_normal")
    t_k = _exact_fraction(flux_tangent, name="flux_tangent")

    actual_normal = magnitude + d_n
    actual_tangent = d_k
    lhs = -(actual_normal * t_n + actual_tangent * t_k)
    leading = -(magnitude * t_n)
    remainder = -(d_n * t_n + d_k * t_k)
    rhs = leading + remainder
    residual = lhs - rhs
    leading_only_defect = lhs - leading

    return ActualShearReplayResult(
        frozen_shear_magnitude=magnitude,
        delta_normal=d_n,
        delta_tangent=d_k,
        flux_normal=t_n,
        flux_tangent=t_k,
        lhs_actual_production=lhs,
        leading_frozen_term=leading,
        retained_remainder=remainder,
        rhs_corrected_split=rhs,
        identity_residual=residual,
        leading_only_defect=leading_only_defect,
    )


def require_corrected_actual_shear_split(**kwargs: ExactInput) -> ActualShearReplayResult:
    """Return the replay result only when the corrected exact split closes."""
    result = replay_actual_shear_identity(**kwargs)
    if not result.exact_identity_verified:
        raise ArithmeticError("corrected actual-shear split did not close exactly")
    return result


def require_leading_only_equivalence(**kwargs: ExactInput) -> ActualShearReplayResult:
    """Fail closed unless omitting the actual-shear remainder is exactly valid.

    This helper makes the corrected-edition distinction executable: a caller
    may use only the frozen leading term precisely when the retained remainder
    is mathematically zero, never because it is numerically small.
    """
    result = require_corrected_actual_shear_split(**kwargs)
    if not result.leading_only_is_valid:
        raise ValueError(
            "leading-only production drops a nonzero actual-shear remainder: "
            f"defect={result.leading_only_defect}"
        )
    return result
