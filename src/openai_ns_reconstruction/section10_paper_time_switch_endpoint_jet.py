"""Endpoint-jet transparency for the pinned Section 10 paper ``timeSwitch``.

The pinned formal source proves that ``timeSwitch`` is identically one for
``t >= 3/4`` and that every positive-order time derivative vanishes for
``t > 3/4``.  Therefore, at ``t = 1``, the finite spacetime jet of
``timeSwitch(t) * F(t,x,y,z)`` is exactly the supplied finite spacetime jet of
``F``: the time-direction Leibniz sum has only its zeroth switch term left.

This module makes that paper-specific analytic reduction executable for scalar
and three-vector jets.  It does *not* construct the missing one-sided Section 9
jet at ``t = 1``, prove that such jets converge to all orders, or extend the
field through ``t = 1``.  The input jet is caller-supplied data; this adapter
only certifies that the official Section 10 time switch is transparent to it at
the endpoint.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .section9_eq921_prefix_jet import MultiIndex4, required_spacetime_multiindices
from .section10_paper_time_switch import Section10PaperTimeSwitchSource


_ENDPOINT_T1 = Fraction(1, 1)


def _nonnegative_order(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("derivative_order must be a nonnegative integer")
    return value


def _required_keys(order: int) -> tuple[MultiIndex4, ...]:
    return required_spacetime_multiindices(order)


def _normalize_scalar_jet(
    values: Mapping[MultiIndex4, object], *, order: int
) -> dict[MultiIndex4, float]:
    if not isinstance(values, Mapping):
        raise TypeError("scalar endpoint jet must be a mapping")
    required = set(_required_keys(order))
    if set(values) != required:
        missing = sorted(required - set(values))
        extra = sorted(set(values) - required)
        raise ValueError(
            "scalar endpoint jet must cover exactly all spacetime derivatives "
            f"through total order {order}; missing={missing}, extra={extra}"
        )
    out: dict[MultiIndex4, float] = {}
    for alpha in _required_keys(order):
        value = values[alpha]
        if isinstance(value, bool):
            raise ValueError(f"scalar endpoint jet[{alpha}] must be finite")
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"scalar endpoint jet[{alpha}] must be finite")
        out[alpha] = number
    return out


def _normalize_vector_jet(
    values: Mapping[MultiIndex4, object], *, order: int
) -> dict[MultiIndex4, np.ndarray]:
    if not isinstance(values, Mapping):
        raise TypeError("vector endpoint jet must be a mapping")
    required = set(_required_keys(order))
    if set(values) != required:
        missing = sorted(required - set(values))
        extra = sorted(set(values) - required)
        raise ValueError(
            "vector endpoint jet must cover exactly all spacetime derivatives "
            f"through total order {order}; missing={missing}, extra={extra}"
        )
    out: dict[MultiIndex4, np.ndarray] = {}
    for alpha in _required_keys(order):
        vector = np.asarray(values[alpha], dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError(f"vector endpoint jet[{alpha}] must be a finite three-vector")
        out[alpha] = np.array(vector, dtype=float, copy=True)
    return out


@dataclass(frozen=True)
class Section10PaperTimeSwitchEndpointJetCertificate:
    """Finite-order endpoint adapter bound to the exact pinned ``timeSwitch``.

    The certificate fixes the endpoint to ``t=1``.  It intentionally accepts no
    arbitrary time window and no caller-supplied switch derivatives.  Those are
    read from the already pinned theorem-facing source: value one and all
    positive time derivatives zero at the endpoint.
    """

    derivative_order: int
    source: Section10PaperTimeSwitchSource

    def __post_init__(self) -> None:
        order = _nonnegative_order(self.derivative_order)
        if not isinstance(self.source, Section10PaperTimeSwitchSource):
            raise TypeError("source must be Section10PaperTimeSwitchSource")
        if self.source != Section10PaperTimeSwitchSource.pinned():
            raise ValueError("source must match the exact pinned Section 10 timeSwitch")
        if not self.source.time_switch_contdiff_theorem_bound:
            raise ValueError("pinned timeSwitch ContDiff theorem is not bound")
        if not self.source.endpoint_t1_switch_value_certified_one:
            raise ValueError("pinned source does not certify timeSwitch(1)=1")
        if not self.source.endpoint_t1_switch_positive_derivatives_certified_zero:
            raise ValueError("pinned source does not certify positive derivatives vanish at t=1")
        object.__setattr__(self, "derivative_order", order)

    @classmethod
    def pinned(cls, derivative_order: int) -> "Section10PaperTimeSwitchEndpointJetCertificate":
        return cls(
            derivative_order=derivative_order,
            source=Section10PaperTimeSwitchSource.pinned(),
        )

    @property
    def endpoint(self) -> Fraction:
        return _ENDPOINT_T1

    def switch_time_derivative(self, order: int) -> Fraction:
        k = _nonnegative_order(order)
        if k > self.derivative_order:
            raise ValueError("requested derivative exceeds certificate order")
        return Fraction(1, 1) if k == 0 else Fraction(0, 1)

    @property
    def switch_time_jet(self) -> tuple[Fraction, ...]:
        return tuple(
            self.switch_time_derivative(k) for k in range(self.derivative_order + 1)
        )

    def apply_scalar_spacetime_jet(
        self, values: Mapping[MultiIndex4, object]
    ) -> Mapping[MultiIndex4, float]:
        """Apply the endpoint Leibniz rule to a supplied scalar spacetime jet."""

        source = _normalize_scalar_jet(values, order=self.derivative_order)
        out: dict[MultiIndex4, float] = {}
        for alpha in _required_keys(self.derivative_order):
            dt, dx, dy, dz = alpha
            total = 0.0
            for k in range(dt + 1):
                beta = (dt - k, dx, dy, dz)
                total += (
                    comb(dt, k)
                    * float(self.switch_time_derivative(k))
                    * source[beta]
                )
            out[alpha] = total
        return MappingProxyType(out)

    def apply_vector_spacetime_jet(
        self, values: Mapping[MultiIndex4, object]
    ) -> Mapping[MultiIndex4, np.ndarray]:
        """Apply the endpoint Leibniz rule to a supplied three-vector spacetime jet."""

        source = _normalize_vector_jet(values, order=self.derivative_order)
        out: dict[MultiIndex4, np.ndarray] = {}
        for alpha in _required_keys(self.derivative_order):
            dt, dx, dy, dz = alpha
            total = np.zeros(3, dtype=float)
            for k in range(dt + 1):
                beta = (dt - k, dx, dy, dz)
                total += (
                    comb(dt, k)
                    * float(self.switch_time_derivative(k))
                    * source[beta]
                )
            total.setflags(write=False)
            out[alpha] = total
        return MappingProxyType(out)

    @property
    def endpoint_multiplier_jet_is_identity(self) -> bool:
        return self.switch_time_jet == (
            Fraction(1, 1),
            *(Fraction(0, 1) for _ in range(self.derivative_order)),
        )

    @property
    def section9_endpoint_jet_supplied_by_actual_construction(self) -> bool:
        return False

    @property
    def section9_field_smooth_extension_through_t1_constructed(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
