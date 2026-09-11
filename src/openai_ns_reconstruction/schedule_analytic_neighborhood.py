"""Constructive complex neighborhood for the actual schedule natural-axis data.

Pinned ``NaturalAxisCoefficients.exists_common_neighborhood`` obtains the
coefficient-space analytic radius by compactness.  This module supplies an
explicit conservative radius and explicit complex-field sup bounds for the
already-landed ``SchedulePressure`` datum.

The construction is deliberately theorem-shaped and fail-closed:

* ``PressureDatum.strip`` has ``|Im z| < 1/2``; we use half that width as a
  deterministic search cap.
* ``complexL`` is kept away from zero by a derivative bound from the real
  interval ``[-11/10,11/10]``.
* ``denominator = H^2 + sigma^2`` is factored as
  ``(H-i sigma)(H+i sigma)``.  On the real interval each factor has modulus at
  least ``sigma``.  A complex-tube perturbation of ``H`` by at most
  ``sigma/4`` therefore leaves both factors at least ``3 sigma/4``.
* The schedule pressure uses ``0 <= shapeExponent <= 1``.  On the tube,
  ``|1+z^2| >= (1-rho)^2`` gives direct bounds for the complex pressure and its
  derivative from an independently certified clock-mass upper bound.
* The eleven actual ``complexField`` inputs then receive explicit sup bounds.
  Their common bound is chosen in the exact shape used by Lean's
  ``finite_family_bound``: ``1 + sum(field_bounds)``.
* Since the tube is convex and contains zero, integrating the complex gradient
  along the straight segment gives a certified upper bound for
  ``sup Re(axisPhase)`` on the same compact tube.

These are executable real-arithmetic inequalities, not interval/Lean proof
objects.  They make the previously opaque ``rho``, common complex-field ``B``,
and ``realPartSup`` inputs explicit for the actual schedule, but do not by
themselves construct the coefficient-space fixed point.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

from .natural_axis import A, D
from .natural_axis_range import NaturalAxisRangeParameters
from .outgoing_tail import TailData
from .schedule_axis_margin import (
    ScheduleLowZMarginWitness,
    certify_schedule_low_Z_margin,
)

_WINDOW_RADIUS = 11.0 / 10.0
_STRIP_HALF_WIDTH = 1.0 / 2.0
_SEARCH_RADIUS = _STRIP_HALF_WIDTH / 2.0
_DENOMINATOR_PERTURBATION_FRACTION = 1.0 / 4.0
_EXPONENT_CAP = 1.0

_FIELD_ORDER = (
    "one",
    "eta",
    "d",
    "inverseL",
    "uStar",
    "uStarEta",
    "wStar",
    "hStar",
    "zStar",
    "chi",
    "gradient",
)


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _finite_nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _h_bound(h: float, j: float, radius: float) -> float:
    """Triangle-inequality bound for ``complexH`` on ``|z| <= radius``."""

    return abs(D(h)) * radius + (1.0 + radius * radius) * (4.0 * radius + abs(j))


def _h_derivative_bound(h: float, j: float, radius: float) -> float:
    """Triangle-inequality bound for the derivative of ``complexH``."""

    return (
        abs(D(h))
        + 2.0 * radius * (4.0 * radius + abs(j))
        + 4.0 * (1.0 + radius * radius)
    )


@dataclass(frozen=True)
class ScheduleAnalyticNeighborhoodCertificate:
    """Explicit tube and field bounds for the actual schedule pressure datum."""

    parameters: NaturalAxisRangeParameters
    sigma: float
    clock_mass_upper: float
    radius: float
    search_radius: float
    complex_abs_radius: float
    l_real_lower: float
    l_tube_lower: float
    h_derivative_upper: float
    h_perturbation_fraction: float
    denominator_lower: float
    pressure_kernel_upper: float
    pressure_upper: float
    pressure_derivative_upper: float
    field_bounds: Mapping[str, float]
    common_field_sup_upper: float
    axis_phase_real_part_sup_upper: float

    def __post_init__(self) -> None:
        sigma = _finite_positive(self.sigma, "sigma")
        mass = _finite_positive(self.clock_mass_upper, "clock_mass_upper")
        radius = _finite_positive(self.radius, "radius")
        search = _finite_positive(self.search_radius, "search_radius")
        if radius > search or search >= _STRIP_HALF_WIDTH:
            raise ValueError("analytic radius must stay strictly inside PressureDatum.strip")
        _finite_positive(self.complex_abs_radius, "complex_abs_radius")
        _finite_positive(self.l_real_lower, "l_real_lower")
        _finite_positive(self.l_tube_lower, "l_tube_lower")
        _finite_positive(self.h_derivative_upper, "h_derivative_upper")
        frac = _finite_nonnegative(
            self.h_perturbation_fraction,
            "h_perturbation_fraction",
        )
        if frac > _DENOMINATOR_PERTURBATION_FRACTION * (1.0 + 1e-12):
            raise ValueError("H perturbation must be at most sigma/4")
        _finite_positive(self.denominator_lower, "denominator_lower")
        _finite_positive(self.pressure_kernel_upper, "pressure_kernel_upper")
        _finite_positive(self.pressure_upper, "pressure_upper")
        _finite_nonnegative(self.pressure_derivative_upper, "pressure_derivative_upper")
        if tuple(self.field_bounds) != _FIELD_ORDER:
            raise ValueError("field_bounds must list the pinned NaturalAxisCoefficients fields")
        for name, value in self.field_bounds.items():
            _finite_positive(value, f"field bound {name}")
        common = _finite_positive(self.common_field_sup_upper, "common_field_sup_upper")
        if common <= max(self.field_bounds.values()):
            raise ValueError("common field bound must strictly dominate every field bound")
        _finite_positive(
            self.axis_phase_real_part_sup_upper,
            "axis_phase_real_part_sup_upper",
        )
        if mass <= 0.0 or sigma <= 0.0:
            raise AssertionError("validated positive inputs must remain positive")

    @property
    def epsilon(self) -> float:
        """Canonical coefficient-space radius ``epsilon = rho/2``."""

        epsilon = self.radius / 2.0
        if epsilon <= 0.0:
            raise ArithmeticError("epsilon underflowed")
        return epsilon


@dataclass(frozen=True)
class ActualScheduleAnalyticInputCertificate:
    """Tie the low-|Z| sigma witness to the explicit analytic neighborhood."""

    low_z: ScheduleLowZMarginWitness
    neighborhood: ScheduleAnalyticNeighborhoodCertificate

    def __post_init__(self) -> None:
        if self.neighborhood.parameters != self.low_z.parameters:
            raise ValueError("low-Z and analytic-neighborhood parameters must agree")
        if self.neighborhood.sigma != self.low_z.cutoff.sigma:
            raise ValueError("analytic neighborhood must use the theorem-selected sigma")
        if self.neighborhood.clock_mass_upper != self.low_z.clock_mass_upper:
            raise ValueError("analytic neighborhood must reuse the certified schedule mass bound")


def certify_axis_analytic_neighborhood(
    *,
    h: float,
    j: float,
    sigma: float,
    clock_mass_upper: float,
) -> ScheduleAnalyticNeighborhoodCertificate:
    """Build explicit ``rho``, complex-field ``B`` and phase-sup certificates.

    ``clock_mass_upper`` must be an independently certified upper bound for
    ``integral clockWeight``.  For the actual schedule this is supplied by
    :func:`schedule_axis_margin.clock_mass_upper`; this function intentionally
    does not estimate it from quadrature samples.
    """

    parameters = NaturalAxisRangeParameters(h=h, j=j)
    sigma = _finite_positive(sigma, "sigma")
    mass = _finite_positive(clock_mass_upper, "clock_mass_upper")

    # Work first on a fixed tube of half the available PressureDatum strip.
    r_search = _SEARCH_RADIUS
    R_search = _WINDOW_RADIUS + r_search

    # L(x)=1-2 h x^2 on the real enlarged window.  NaturalAxisRange guarantees
    # h>0 and h<=1/100, so this lower bound is positive.
    l_real_lower = 1.0 - 2.0 * parameters.h * _WINDOW_RADIUS**2
    l_real_lower = _finite_positive(l_real_lower, "real-axis L lower bound")
    l_derivative_upper = 4.0 * parameters.h * R_search
    l_derivative_upper = _finite_positive(
        l_derivative_upper,
        "complexL derivative upper bound",
    )
    # Leave at least 3/4 of the real-axis L gap.
    r_l = l_real_lower / (4.0 * l_derivative_upper)

    # If z is within rho of a real x, then
    # |H(z)-H(x)| <= sup|H'| rho.  Restrict this to sigma/4.  Since H(x) is
    # real, |H(x) +/- i sigma| >= sigma and both complex factors remain >=3sigma/4.
    h_derivative_upper = _finite_positive(
        _h_derivative_bound(parameters.h, parameters.j, R_search),
        "complexH derivative upper bound",
    )
    r_den = (
        _DENOMINATOR_PERTURBATION_FRACTION
        * sigma
        / h_derivative_upper
    )
    r_den = _finite_positive(r_den, "denominator-safe radius")

    radius = min(r_search, r_l, r_den)
    radius = _finite_positive(radius, "common analytic radius")
    R = _WINDOW_RADIUS + radius

    l_tube_lower = l_real_lower - l_derivative_upper * radius
    l_tube_lower = _finite_positive(l_tube_lower, "tube L lower bound")

    h_perturbation_fraction = h_derivative_upper * radius / sigma
    if h_perturbation_fraction > _DENOMINATOR_PERTURBATION_FRACTION * (1.0 + 1e-12):
        raise ArithmeticError("floating realization lost the H perturbation certificate")
    factor_lower = sigma * (1.0 - h_perturbation_fraction)
    denominator_lower = _finite_positive(
        factor_lower * factor_lower,
        "denominator lower bound",
    )

    # On |Im z|<=rho, 1+z^2=(z-i)(z+i), so both factors are at least 1-rho.
    base_lower = (1.0 - radius) ** 2
    base_lower = _finite_positive(base_lower, "|1+z^2| lower bound")
    pressure_kernel_upper = _finite_positive(
        base_lower ** (-2.0 * _EXPONENT_CAP),
        "complex pressure-kernel upper bound",
    )
    pressure_upper = _finite_positive(
        0.5 * mass * pressure_kernel_upper,
        "complex pressure upper bound",
    )
    pressure_derivative_upper = _finite_nonnegative(
        2.0
        * mass
        * _EXPONENT_CAP
        * R
        * pressure_kernel_upper
        / base_lower,
        "complex pressure-derivative upper bound",
    )

    d_bound = 1.0 + R * R
    l_abs_bound = 1.0 + 2.0 * parameters.h * R * R
    u_bound = 4.0 * R + abs(parameters.j)
    h_bound = _h_bound(parameters.h, parameters.j, R)
    w_bound = (
        1.0
        + 4.0 * d_bound
        + 2.0 * abs(D(parameters.h)) * R * u_bound
    )

    # Use the factor separation instead of the much looser global
    # numerator/denominator quotient.  With c=|H-H(real)|/sigma<1,
    # |H^2/(H^2+sigma^2)| <= 1 + (1-c)^-2 and
    # |H/(H^2+sigma^2)| <= 1/(sigma(1-c)).
    factor_fraction = 1.0 - h_perturbation_fraction
    chi_bound = 1.0 + 1.0 / (factor_fraction * factor_fraction)
    gradient_bound = l_abs_bound / (sigma * factor_fraction)

    a = abs(A(parameters.h))
    z_bound = (
        a * (1.0 + 2.0 * R * u_bound) * u_bound
        + 4.0 * h_bound
        + d_bound * pressure_derivative_upper
        + 4.0 * a * R * pressure_upper
    )

    bounds_dict = {
        "one": 1.0,
        "eta": R,
        "d": d_bound,
        "inverseL": 1.0 / l_tube_lower,
        "uStar": u_bound,
        "uStarEta": 4.0,
        "wStar": w_bound,
        "hStar": h_bound,
        "zStar": z_bound,
        "chi": chi_bound,
        "gradient": gradient_bound,
    }
    bounds = MappingProxyType(
        {
            name: _finite_positive(bounds_dict[name], f"field bound {name}")
            for name in _FIELD_ORDER
        }
    )

    # Match the constructive positive common bound in Lean finite_family_bound.
    common = 1.0 + math.fsum(bounds.values())
    common = _finite_positive(common, "common complex-field sup bound")

    # The closed tube is convex and contains 0.  The straight segment from 0 to
    # z stays inside it, so |axisPhase(z)| <= |z| sup|complexGradient|.
    phase_sup = _finite_positive(
        R * bounds["gradient"],
        "axisPhase real-part sup upper bound",
    )

    return ScheduleAnalyticNeighborhoodCertificate(
        parameters=parameters,
        sigma=sigma,
        clock_mass_upper=mass,
        radius=radius,
        search_radius=r_search,
        complex_abs_radius=R,
        l_real_lower=l_real_lower,
        l_tube_lower=l_tube_lower,
        h_derivative_upper=h_derivative_upper,
        h_perturbation_fraction=h_perturbation_fraction,
        denominator_lower=denominator_lower,
        pressure_kernel_upper=pressure_kernel_upper,
        pressure_upper=pressure_upper,
        pressure_derivative_upper=pressure_derivative_upper,
        field_bounds=bounds,
        common_field_sup_upper=common,
        axis_phase_real_part_sup_upper=phase_sup,
    )


def certify_actual_schedule_analytic_inputs(
    data: TailData,
    j: float,
) -> ActualScheduleAnalyticInputCertificate:
    """Instantiate the analytic-neighborhood certificate for ``SchedulePressure``.

    The theorem-selected ``sigma`` and the no-sampling clock-mass upper bound are
    both reused from the existing low-|Z| schedule witness.  No free numerical
    radius or complex sup bound is supplied by the caller.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    low_z = certify_schedule_low_Z_margin(data, j)
    neighborhood = certify_axis_analytic_neighborhood(
        h=data.h,
        j=j,
        sigma=low_z.cutoff.sigma,
        clock_mass_upper=low_z.clock_mass_upper,
    )
    return ActualScheduleAnalyticInputCertificate(
        low_z=low_z,
        neighborhood=neighborhood,
    )
