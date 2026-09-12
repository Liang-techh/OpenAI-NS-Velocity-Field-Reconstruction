"""Pinned coefficient product ``u1 * u1`` at the genuine first Picard state.

The landed first Picard axial field is represented without collapsing theorem
scales as

    u1 = r + b / Lambda + c / Lambda^2 + a^2 p / Lambda,

coefficient-jet by coefficient-jet.  Here ``r,b,c`` are the ordinary Decimal
pieces already carried by ``ActualScheduleWideFirstPicardState`` and ``p`` is
one half of the landed wide-pressure normalized factor.  The common ``a^2``
scale is the genuine theorem-selected

    a(eta) = exp(Lambda * realPhase(eta)) / C,

not a surrogate amplitude.

This module lifts exactly one nonlinear ``AxisOperators.product`` seam needed by
``naturalRemainder(x1)``: ``u1 * u1``.  The pinned radial convolution and
eta-Leibniz rule are applied before collecting scales.  The result is kept as

* ordinary powers Lambda^0 through Lambda^-4;
* pressure-linear ``a^2`` powers Lambda^-1 through Lambda^-3; and
* the pressure-square ``a^4`` power Lambda^-2.

No O(1) value is ever added to a tiny inverse-Lambda correction, and the
``a^2``/``a^4`` pieces are exposed as signed-log jets only after their Decimal
normalized numerators are formed.  No caller-supplied Lambda, amplitude,
pressure table, coefficient table, derivative table, or product cutoff is
accepted.

This is one signed-log-aware product primitive, not the complete
``naturalRemainder(x1)``, ``x2``, a fixed-point convergence certificate, global
weighted ``AxisSpace`` membership, or a paper-exact velocity.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_MAX_ORDINARY_POWER = 4
_TWO = Decimal(2)


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not WINDOW_LEFT <= eta <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _signed_log_term(
    numerator: Decimal,
    *,
    amplitude_log: Decimal,
    amplitude_power: int,
    Lambda: Decimal,
    inverse_lambda_power: int,
) -> SignedLogCoefficientJet:
    """Encode ``a^amplitude_power * numerator / Lambda^inverse_lambda_power``.

    The enormous common amplitude logarithm is kept in ``log_scale`` while the
    moderate numerator/Lambda correction stays in ``log_factor``.  This mirrors
    the existing signed-log amplitude/pressure representation and avoids adding
    quantities whose decimal exponents differ by hundreds of orders.
    """

    _finite_decimal(numerator, "numerator")
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if amplitude_power <= 0 or inverse_lambda_power < 0:
        raise ValueError("invalid signed-log scale powers")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()

    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_scale = +(Decimal(amplitude_power) * amplitude_log)
        log_factor = +(
            abs(numerator).ln()
            - Decimal(inverse_lambda_power) * Lambda.ln()
        )
    return SignedLogCoefficientJet(
        sign=1 if numerator > 0 else -1,
        log_scale=log_scale,
        log_factor=log_factor,
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardAxialSquareCoefficientJet:
    """One coefficient jet of ``u1 * u1`` with all theorem scales split.

    The represented coefficient is

    ``o0 + o1/L + o2/L^2 + o3/L^3 + o4/L^4``
    ``+ a^2 (q1/L + q2/L^2 + q3/L^3)``
    ``+ a^4 s2/L^2``.

    ``a`` is the actual theorem-selected amplitude at the requested ``eta``.
    The class deliberately exposes no collapsed total-value method.
    """

    ordinary_reference: Decimal
    ordinary_inverse_lambda_numerator: Decimal
    ordinary_inverse_lambda_squared_numerator: Decimal
    ordinary_inverse_lambda_cubed_numerator: Decimal
    ordinary_inverse_lambda_fourth_numerator: Decimal
    pressure_linear_inverse_lambda_numerator: Decimal
    pressure_linear_inverse_lambda_squared_numerator: Decimal
    pressure_linear_inverse_lambda_cubed_numerator: Decimal
    pressure_square_inverse_lambda_squared_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in (
            ("ordinary_reference", self.ordinary_reference),
            ("ordinary_inverse_lambda_numerator", self.ordinary_inverse_lambda_numerator),
            (
                "ordinary_inverse_lambda_squared_numerator",
                self.ordinary_inverse_lambda_squared_numerator,
            ),
            (
                "ordinary_inverse_lambda_cubed_numerator",
                self.ordinary_inverse_lambda_cubed_numerator,
            ),
            (
                "ordinary_inverse_lambda_fourth_numerator",
                self.ordinary_inverse_lambda_fourth_numerator,
            ),
            (
                "pressure_linear_inverse_lambda_numerator",
                self.pressure_linear_inverse_lambda_numerator,
            ),
            (
                "pressure_linear_inverse_lambda_squared_numerator",
                self.pressure_linear_inverse_lambda_squared_numerator,
            ),
            (
                "pressure_linear_inverse_lambda_cubed_numerator",
                self.pressure_linear_inverse_lambda_cubed_numerator,
            ),
            (
                "pressure_square_inverse_lambda_squared_numerator",
                self.pressure_square_inverse_lambda_squared_numerator,
            ),
            ("Lambda", self.Lambda),
            ("amplitude_log", self.amplitude_log),
        ):
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_correction_terms_decimal(
        self,
    ) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        """Return ordinary Lambda^-1 through Lambda^-4 terms separately."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            inverse2 = +(inverse * inverse)
            return (
                +(self.ordinary_inverse_lambda_numerator * inverse),
                +(self.ordinary_inverse_lambda_squared_numerator * inverse2),
                +(self.ordinary_inverse_lambda_cubed_numerator * inverse2 * inverse),
                +(self.ordinary_inverse_lambda_fourth_numerator * inverse2 * inverse2),
            )

    def pressure_linear_terms_log(
        self,
    ) -> tuple[SignedLogCoefficientJet, SignedLogCoefficientJet, SignedLogCoefficientJet]:
        """Return the ``a^2`` pressure-linear Lambda^-1..-3 terms in signed-log form."""

        return (
            _signed_log_term(
                self.pressure_linear_inverse_lambda_numerator,
                amplitude_log=self.amplitude_log,
                amplitude_power=2,
                Lambda=self.Lambda,
                inverse_lambda_power=1,
            ),
            _signed_log_term(
                self.pressure_linear_inverse_lambda_squared_numerator,
                amplitude_log=self.amplitude_log,
                amplitude_power=2,
                Lambda=self.Lambda,
                inverse_lambda_power=2,
            ),
            _signed_log_term(
                self.pressure_linear_inverse_lambda_cubed_numerator,
                amplitude_log=self.amplitude_log,
                amplitude_power=2,
                Lambda=self.Lambda,
                inverse_lambda_power=3,
            ),
        )

    def pressure_square_term_log(self) -> SignedLogCoefficientJet:
        """Return the ``a^4 / Lambda^2`` pressure-square term in signed-log form."""

        return _signed_log_term(
            self.pressure_square_inverse_lambda_squared_numerator,
            amplitude_log=self.amplitude_log,
            amplitude_power=4,
            Lambda=self.Lambda,
            inverse_lambda_power=2,
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardAxialSquareState:
    """The pinned product ``u1 * u1`` on actual theorem-selected schedule data."""

    x1: ActualScheduleWideFirstPicardState

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")
        if self.x1.Lambda != self.x1.remainder.axial.wide_pressure.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual theorem-selected amplitude scale")

    @property
    def Lambda(self) -> Decimal:
        return self.x1.Lambda

    @property
    def epsilon(self) -> float:
        return self.x1.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def mixed_scale_axial_square_materialized(self) -> bool:
        return True

    @property
    def signed_log_pressure_product_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    def _pressure_inverse_lambda_numerator(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Return ``p`` where the x1 pressure jet is ``a^2 * p / Lambda``.

        The landed axial remainder exposes ``q`` such that its pressure jet is
        ``a^2 q``.  The outer Picard scale contributes ``1/(2 Lambda)``, hence
        ``p=q/2``.  Keeping ``Lambda`` outside is what lets subsequent products
        preserve every inverse-Lambda power exactly.
        """

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(
                self.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
                / _TWO
            )

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardAxialSquareCoefficientJet:
        """Evaluate the exact pinned radial-convolution/eta-Leibniz square jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        # Force the landed x1 guards at the requested coordinate before the
        # convolution.  All summands below remain tied to this same state.
        self.x1.jet_pair(n, m, eta)

        ordinary = [Decimal(0) for _ in range(_MAX_ORDINARY_POWER + 1)]
        pressure_linear = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
        pressure_square_power2 = Decimal(0)

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    l = m - k
                    _, left = self.x1.jet_pair(i, k, eta)
                    _, right = self.x1.jet_pair(j, l, eta)
                    weight = Decimal(math.comb(m, k))

                    left_ordinary = (
                        left.reference,
                        left.inverse_lambda_numerator,
                        left.inverse_lambda_squared_numerator,
                    )
                    right_ordinary = (
                        right.reference,
                        right.inverse_lambda_numerator,
                        right.inverse_lambda_squared_numerator,
                    )
                    for left_power, left_value in enumerate(left_ordinary):
                        for right_power, right_value in enumerate(right_ordinary):
                            ordinary[left_power + right_power] += (
                                weight * left_value * right_value
                            )

                    left_pressure = self._pressure_inverse_lambda_numerator(
                        i, k, eta
                    )
                    right_pressure = self._pressure_inverse_lambda_numerator(
                        j, l, eta
                    )
                    for left_power, left_value in enumerate(left_ordinary):
                        pressure_linear[left_power + 1] += (
                            weight * left_value * right_pressure
                        )
                    for right_power, right_value in enumerate(right_ordinary):
                        pressure_linear[right_power + 1] += (
                            weight * left_pressure * right_value
                        )
                    pressure_square_power2 += (
                        weight * left_pressure * right_pressure
                    )

            ordinary = [+value for value in ordinary]
            pressure_linear = [+value for value in pressure_linear]
            pressure_square_power2 = +pressure_square_power2

        amplitude_log = (
            self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
        )
        return MixedScaleFirstPicardAxialSquareCoefficientJet(
            ordinary_reference=ordinary[0],
            ordinary_inverse_lambda_numerator=ordinary[1],
            ordinary_inverse_lambda_squared_numerator=ordinary[2],
            ordinary_inverse_lambda_cubed_numerator=ordinary[3],
            ordinary_inverse_lambda_fourth_numerator=ordinary[4],
            pressure_linear_inverse_lambda_numerator=pressure_linear[1],
            pressure_linear_inverse_lambda_squared_numerator=pressure_linear[2],
            pressure_linear_inverse_lambda_cubed_numerator=pressure_linear[3],
            pressure_square_inverse_lambda_squared_numerator=pressure_square_power2,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_axial_square_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardAxialSquareState:
    """Lift the pinned product to ``u1 * u1`` without theorem-scale collapse."""

    return ActualScheduleWideFirstPicardAxialSquareState(x1=x1)
