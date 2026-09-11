"""Conditional endpoint-limit certificate for Section 10 residual jets.

``CandidateFromLimits.lean`` assumes that every full spacetime derivative of
an actual closed-past Navier--Stokes residual has a locally uniform limit as
``t -> 1-``.  The repository already has adapters for *using* those endpoint
jets, but until the upstream field is complete it cannot honestly manufacture
the limits themselves.

This module records a reusable sufficient condition for that missing step.  On
one compact spatial window and for one full derivative degree, suppose an
independently proved estimate gives

    sup_x || d/dt D^n R(t,x) || <= C (T-t)^(-alpha),

for ``valid_from <= t < T``, with ``C >= 0`` and ``0 <= alpha < 1``.  Then the
fundamental theorem of calculus gives the uniform Cauchy modulus

    sup_x ||D^n R(s,x)-D^n R(t,x)||
      <= C/(1-alpha) * ((T-t)^(1-alpha) - (T-s)^(1-alpha)),

and the endpoint tail bound

    sup_x ||L_n(x)-D^n R(t,x)||
      <= C/(1-alpha) * (T-t)^(1-alpha).

Only this implication is implemented here.  The majorant itself is caller
supplied and must eventually be derived from the actual localized residual;
finite sampling, fitted slopes, or a manufactured field do not satisfy that
hypothesis.  Consequently this is formal-structure infrastructure and cannot
promote ``paper_exact_velocity_available``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import operator


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(out)


def _finite_real(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real number")
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} must be a finite real number")
    return out


@dataclass(frozen=True)
class UniformCauchyCertificate:
    """Arithmetic consequence of one supplied endpoint derivative majorant.

    ``certified`` means only that the recorded numbers obey the integrated
    majorant arithmetic.  It does *not* certify that the actual residual obeys
    the majorant hypothesis.
    """

    derivative_degree: int
    spatial_window: int
    t0: float
    t1: float
    endpoint: float
    interval_upper_bound: float
    endpoint_tail_upper_bound: float

    @property
    def certified(self) -> bool:
        values = (
            self.t0,
            self.t1,
            self.endpoint,
            self.interval_upper_bound,
            self.endpoint_tail_upper_bound,
        )
        return (
            self.derivative_degree >= 0
            and self.spatial_window >= 0
            and all(math.isfinite(value) for value in values)
            and self.t0 <= self.t1 < self.endpoint
            and 0.0 <= self.interval_upper_bound
            <= self.endpoint_tail_upper_bound
        )


@dataclass(frozen=True)
class EndpointPowerLawMajorant:
    """Supplied integrable power-law bound for one residual-jet derivative.

    Semantically the caller is asserting the *independently proved* bound

    ``sup_{x in K_m} ||partial_t D^n R(t,x)|| <= C (T-t)^(-alpha)``

    on ``valid_from <= t < T``.  This class validates the parameters and
    performs only the exact antiderivative arithmetic implied by that bound.
    """

    derivative_degree: int
    spatial_window: int
    coefficient: float
    singularity_exponent: float
    valid_from: float = 0.75
    endpoint: float = 1.0

    def __post_init__(self) -> None:
        degree = _natural(self.derivative_degree, "derivative degree")
        window = _natural(self.spatial_window, "spatial window")
        coefficient = _finite_real(self.coefficient, "coefficient")
        exponent = _finite_real(self.singularity_exponent, "singularity exponent")
        valid_from = _finite_real(self.valid_from, "valid_from")
        endpoint = _finite_real(self.endpoint, "endpoint")
        if coefficient < 0.0:
            raise ValueError("coefficient must be nonnegative")
        if not 0.0 <= exponent < 1.0:
            raise ValueError("singularity exponent must satisfy 0 <= alpha < 1")
        if valid_from >= endpoint:
            raise ValueError("valid_from must be strictly before endpoint")
        object.__setattr__(self, "derivative_degree", degree)
        object.__setattr__(self, "spatial_window", window)
        object.__setattr__(self, "coefficient", coefficient)
        object.__setattr__(self, "singularity_exponent", exponent)
        object.__setattr__(self, "valid_from", valid_from)
        object.__setattr__(self, "endpoint", endpoint)

    @property
    def integrability_power(self) -> float:
        return 1.0 - self.singularity_exponent

    def _time(self, value: object, name: str) -> float:
        time = _finite_real(value, name)
        if time < self.valid_from or time >= self.endpoint:
            raise ValueError(
                f"{name} must satisfy valid_from <= {name} < endpoint"
            )
        return time

    def derivative_upper_bound(self, t: object) -> float:
        """Return ``C (T-t)^(-alpha)`` for a time in the certified regime."""
        t = self._time(t, "t")
        if self.coefficient == 0.0:
            return 0.0
        value = self.coefficient * math.pow(
            self.endpoint - t, -self.singularity_exponent
        )
        if not math.isfinite(value):
            raise ArithmeticError("endpoint derivative majorant overflowed")
        return value

    def endpoint_tail_bound(self, t: object) -> float:
        """Uniform distance-to-limit budget from time ``t`` to the endpoint."""
        t = self._time(t, "t")
        if self.coefficient == 0.0:
            return 0.0
        power = self.integrability_power
        value = (
            self.coefficient
            / power
            * math.pow(self.endpoint - t, power)
        )
        if not math.isfinite(value):
            raise ArithmeticError("endpoint tail majorant overflowed")
        return value

    def interval_bound(self, t0: object, t1: object) -> float:
        """Uniform Cauchy budget between two pre-endpoint times."""
        t0 = self._time(t0, "t0")
        t1 = self._time(t1, "t1")
        if t1 < t0:
            raise ValueError("interval bound requires t0 <= t1")
        if self.coefficient == 0.0 or t0 == t1:
            return 0.0
        power = self.integrability_power
        value = self.coefficient / power * (
            math.pow(self.endpoint - t0, power)
            - math.pow(self.endpoint - t1, power)
        )
        if not math.isfinite(value):
            raise ArithmeticError("endpoint interval majorant overflowed")
        if value < 0.0:
            raise ArithmeticError("endpoint interval majorant became negative")
        return value

    def cauchy_certificate(self, t0: object, t1: object) -> UniformCauchyCertificate:
        """Package one conditional locally-uniform Cauchy implication."""
        t0_value = self._time(t0, "t0")
        t1_value = self._time(t1, "t1")
        if t1_value < t0_value:
            raise ValueError("Cauchy certificate requires t0 <= t1")
        interval = self.interval_bound(t0_value, t1_value)
        tail = self.endpoint_tail_bound(t0_value)
        if interval > tail and not math.isclose(interval, tail, rel_tol=8e-16, abs_tol=0.0):
            raise ArithmeticError("interval majorant exceeds its endpoint tail budget")
        certificate = UniformCauchyCertificate(
            derivative_degree=self.derivative_degree,
            spatial_window=self.spatial_window,
            t0=t0_value,
            t1=t1_value,
            endpoint=self.endpoint,
            interval_upper_bound=interval,
            endpoint_tail_upper_bound=tail,
        )
        if not certificate.certified:
            raise ArithmeticError("endpoint Cauchy certificate invariant failed")
        return certificate

    def epsilon_entry_time(self, epsilon: object) -> float:
        """Earliest time after which the endpoint tail budget is <= epsilon.

        This is a modulus-of-convergence calculation, conditional on the same
        independently certified derivative majorant.  It never estimates the
        coefficient or exponent from residual samples.
        """
        epsilon = _finite_real(epsilon, "epsilon")
        if epsilon <= 0.0:
            raise ValueError("epsilon must be positive")
        if self.coefficient == 0.0:
            return self.valid_from
        if self.endpoint_tail_bound(self.valid_from) <= epsilon:
            return self.valid_from
        power = self.integrability_power
        radius = math.pow(epsilon * power / self.coefficient, 1.0 / power)
        if not math.isfinite(radius) or radius <= 0.0:
            raise ArithmeticError("epsilon endpoint window is not representable")
        entry = self.endpoint - radius
        if entry >= self.endpoint:
            raise ArithmeticError("epsilon entry time is below floating resolution")
        return max(self.valid_from, entry)
