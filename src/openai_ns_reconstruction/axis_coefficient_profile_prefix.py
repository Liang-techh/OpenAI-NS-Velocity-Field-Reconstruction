"""Finite radial polynomial prefixes of the formal mixed-scale coefficient state.

The formal solver supplies coefficient rows, not a converged profile.  This
module evaluates a finite natural-radial polynomial prefix while preserving
every sparse ``(amplitude_power, inverse_Lambda_power)`` channel.  It never
forms an amplitude or Lambda power and it does not infer a tail bound from the
retained rows.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import math
from types import MappingProxyType

from .axis_amplitude_log_scale import AmplitudeLogSource
from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_mixed_scale import (
    Channel,
    MixedScaleCoefficient,
)
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_formal_solver import (
    FormalAxisCoefficientSolverState,
)


_DECIMAL_PRECISION = 96


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not WINDOW_LEFT <= eta <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _decimal_y(value: Decimal | int) -> Decimal:
    if isinstance(value, bool):
        raise TypeError("Y must be a finite nonnegative Decimal or integer")
    if isinstance(value, int):
        value = Decimal(value)
    if not isinstance(value, Decimal) or not value.is_finite() or value < 0:
        raise ValueError("Y must be a finite nonnegative Decimal")
    return value


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _falling_factorial(n: int, order: int) -> Decimal:
    if order > n:
        return Decimal(0)
    return Decimal(math.factorial(n) // math.factorial(n - order))


def evaluate_mixed_scale_radial_prefix(
    coefficients: Sequence[MixedScaleCoefficient],
    Y: Decimal | int,
    radial_order: int = 0,
    *,
    average: bool = False,
) -> MixedScaleCoefficient:
    """Evaluate one finite radial derivative/average map on coefficient rows.

    For row ``n`` the contribution is

    ``n!/(n-r)! * Y^(n-r) * coefficients[n]``

    when ``n >= r``.  With ``average=True`` the original coefficient is first
    divided by ``n + 1``.  The helper accepts finite mathematical prefixes and
    does not identify them with an actual theorem field.
    """

    radial_order = _index(radial_order, "radial_order")
    Y = _decimal_y(Y)
    if not isinstance(coefficients, Sequence):
        raise TypeError("coefficients must be a finite sequence of channel maps")
    if not isinstance(average, bool):
        raise TypeError("average must be a boolean")
    rows = tuple(coefficients)
    for row in rows:
        if not isinstance(row, MixedScaleCoefficient):
            raise TypeError("each coefficient row must be MixedScaleCoefficient")

    total = MixedScaleCoefficient.zero()
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        for n, coefficient in enumerate(rows):
            if n < radial_order:
                continue
            exponent = n - radial_order
            if Y == 0 and exponent > 0:
                continue
            y_power = Decimal(1) if exponent == 0 else Y**exponent
            scalar = _falling_factorial(n, radial_order) * y_power
            source_coefficient = coefficient
            if average:
                source_coefficient = coefficient.scale(+(
                    Decimal(1) / Decimal(n + 1)
                ))
            total = total.add(source_coefficient.scale(+scalar))
    return total


def _signed_log_terms(
    coefficient: MixedScaleCoefficient,
    *,
    amplitude_log: Decimal,
    Lambda: Decimal,
    amplitude_log_source: AmplitudeLogSource | None = None,
) -> Mapping[Channel, SignedLogCoefficientJet]:
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if amplitude_log_source is not None and not isinstance(
        amplitude_log_source,
        AmplitudeLogSource,
    ):
        raise TypeError("amplitude_log_source must be AmplitudeLogSource or None")
    if amplitude_log_source is not None:
        if amplitude_log_source.midpoint != amplitude_log:
            raise ValueError(
                "amplitude_log_source midpoint must match amplitude_log"
            )
        if amplitude_log_source.enclosure.Lambda != Lambda:
            raise ValueError("amplitude_log_source Lambda must match Lambda")
    values: dict[Channel, SignedLogCoefficientJet] = {}
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        for (q, p), numerator in coefficient.items():
            if numerator == 0:
                continue
            # Keep the finite Decimal amplitude log intact when multiplying by
            # a channel power. The factor log remains at normal local
            # precision; the two components stay separate in the result.
            if q == 0:
                log_scale = Decimal(0)
            else:
                scale_precision = max(
                    _DECIMAL_PRECISION,
                    len(amplitude_log.as_tuple().digits) + len(str(q)),
                )
                with localcontext() as scale_ctx:
                    scale_ctx.prec = scale_precision
                    log_scale = +(Decimal(q) * amplitude_log)
            amplitude_log_scale = (
                None
                if amplitude_log_source is None
                else amplitude_log_source.power(q, log_scale)
            )
            values[(q, p)] = SignedLogCoefficientJet(
                sign=1 if numerator > 0 else -1,
                log_scale=log_scale,
                log_factor=+(
                    abs(numerator).ln() - Decimal(p) * Lambda.ln()
                ),
                amplitude_log_scale=amplitude_log_scale,
            )
    return MappingProxyType(values)


@dataclass(frozen=True)
class FormalAxisProfilePrefix:
    """One finite formal radial prefix for angular, axial, average, and pressure maps."""

    angular: MixedScaleCoefficient
    axial: MixedScaleCoefficient
    axial_average: MixedScaleCoefficient
    pressure: MixedScaleCoefficient
    Lambda: Decimal
    amplitude_log: Decimal
    max_n: int
    radial_order: int
    eta_order: int
    Y: Decimal
    eta: float
    amplitude_log_source: AmplitudeLogSource | None = None

    def __post_init__(self) -> None:
        for name in ("angular", "axial", "axial_average", "pressure"):
            if not isinstance(getattr(self, name), MixedScaleCoefficient):
                raise TypeError(f"{name} must be MixedScaleCoefficient")
        _finite_decimal(self.Lambda, "Lambda")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")
        _finite_decimal(self.amplitude_log, "amplitude_log")
        _index(self.max_n, "max_n")
        _index(self.radial_order, "radial_order")
        _index(self.eta_order, "eta_order")
        _decimal_y(self.Y)
        eta = _eta_in_window(self.eta)
        source = self.amplitude_log_source
        if source is not None:
            if not isinstance(source, AmplitudeLogSource):
                raise TypeError(
                    "amplitude_log_source must be AmplitudeLogSource or None"
                )
            if source.midpoint != self.amplitude_log:
                raise ValueError(
                    "amplitude_log_source midpoint must match amplitude_log"
                )
            if source.enclosure.Lambda != self.Lambda:
                raise ValueError(
                    "amplitude_log_source Lambda must match Lambda"
                )
            if source.eta != Fraction.from_float(eta):
                raise ValueError(
                    "amplitude_log_source eta must match the profile eta"
                )

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def truncation_certified(self) -> bool:
        return False

    def angular_terms_log(self) -> Mapping[Channel, SignedLogCoefficientJet]:
        """Return signed-log views of the angular sparse channels."""

        return _signed_log_terms(
            self.angular,
            amplitude_log=self.amplitude_log,
            Lambda=self.Lambda,
            amplitude_log_source=self.amplitude_log_source,
        )

    def axial_terms_log(self) -> Mapping[Channel, SignedLogCoefficientJet]:
        """Return signed-log views of the axial sparse channels."""

        return _signed_log_terms(
            self.axial,
            amplitude_log=self.amplitude_log,
            Lambda=self.Lambda,
            amplitude_log_source=self.amplitude_log_source,
        )

    def axial_average_terms_log(self) -> Mapping[Channel, SignedLogCoefficientJet]:
        """Return signed-log views of the averaged axial sparse channels."""

        return _signed_log_terms(
            self.axial_average,
            amplitude_log=self.amplitude_log,
            Lambda=self.Lambda,
            amplitude_log_source=self.amplitude_log_source,
        )

    def pressure_terms_log(self) -> Mapping[Channel, SignedLogCoefficientJet]:
        """Return signed-log views of the scaled pressure primitive channels."""

        return _signed_log_terms(
            self.pressure,
            amplitude_log=self.amplitude_log,
            Lambda=self.Lambda,
            amplitude_log_source=self.amplitude_log_source,
        )


def _actual_amplitude_log_source(
    solver: FormalAxisCoefficientSolverState,
    eta: float,
) -> AmplitudeLogSource:
    """Return the source bound to the actual x1 amplitude, fail-closed."""

    try:
        amplitude = solver.x1.remainder.axial.wide_pressure.amplitude
    except AttributeError as error:
        raise ValueError(
            "actual solver is missing its bound amplitude log source"
        ) from error
    source_factory = getattr(amplitude, "log_amplitude_source", None)
    if not callable(source_factory):
        raise ValueError(
            "actual amplitude does not expose log_amplitude_source"
        )
    source = source_factory(eta)
    if not isinstance(source, AmplitudeLogSource):
        raise TypeError("actual amplitude log source has the wrong type")
    if source.eta != Fraction.from_float(eta):
        raise ValueError("actual amplitude log source eta mismatch")
    return source


def formal_axis_profile_prefix(
    solver: FormalAxisCoefficientSolverState,
    max_n: int,
    Y: Decimal | int,
    eta: float,
    radial_order: int = 0,
    eta_order: int = 0,
) -> FormalAxisProfilePrefix:
    """Evaluate a finite formal radial prefix anchored by ``solver``."""

    if not isinstance(solver, FormalAxisCoefficientSolverState):
        raise TypeError("solver must be FormalAxisCoefficientSolverState")
    max_n = _index(max_n, "max_n")
    radial_order = _index(radial_order, "radial_order")
    eta_order = _index(eta_order, "eta_order")
    Y = _decimal_y(Y)
    eta = _eta_in_window(eta)
    amplitude_log_source = _actual_amplitude_log_source(solver, eta)
    rows = solver.profile_jet_prefix(max_n, eta_order, eta)
    angular_rows = tuple(pair[0] for pair in rows)
    axial_rows = tuple(pair[1] for pair in rows)
    pressure_rows = tuple(pair[2] for pair in rows)
    angular = evaluate_mixed_scale_radial_prefix(
        angular_rows,
        Y,
        radial_order,
    )
    axial = evaluate_mixed_scale_radial_prefix(
        axial_rows,
        Y,
        radial_order,
    )
    axial_average = evaluate_mixed_scale_radial_prefix(
        axial_rows,
        Y,
        radial_order,
        average=True,
    )
    pressure = evaluate_mixed_scale_radial_prefix(
        pressure_rows,
        Y,
        radial_order,
    )
    return FormalAxisProfilePrefix(
        angular=angular,
        axial=axial,
        axial_average=axial_average,
        pressure=pressure,
        Lambda=solver.Lambda,
        amplitude_log=amplitude_log_source.midpoint,
        max_n=max_n,
        radial_order=radial_order,
        eta_order=eta_order,
        Y=Y,
        eta=eta,
        amplitude_log_source=amplitude_log_source,
    )


__all__ = [
    "FormalAxisProfilePrefix",
    "evaluate_mixed_scale_radial_prefix",
    "formal_axis_profile_prefix",
]
