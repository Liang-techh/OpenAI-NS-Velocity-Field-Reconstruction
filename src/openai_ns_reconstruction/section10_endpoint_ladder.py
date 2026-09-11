"""Finite-order Section 10 endpoint-majorant ladder.

The pinned Section 10 endpoint construction needs locally uniform limits for
*every* spacetime derivative of the closed-past residual as ``t -> 1-``.
``endpoint_limit_majorant.py`` records a sufficient integrable time-derivative
majorant for one derivative degree, while ``time_localization.py`` fail-closes
that majorant to the official ``t >= 3/4`` unit plateau.

This module packages a contiguous finite family of those hypotheses.  It is a
finite-order readiness certificate only: it does not derive the majorants from
the actual residual, manufacture endpoint jets, or prove all-order smoothness.
In particular, it must never promote ``paper_exact_velocity_available``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

from .endpoint_limit_majorant import EndpointPowerLawMajorant, UniformCauchyCertificate
from .time_localization import (
    LATE_START,
    SECTION10_ENDPOINT,
    EndpointLocalizationTransferCertificate,
    section10_endpoint_localization_transfer,
)


@dataclass(frozen=True)
class Section10EndpointMajorantLadder:
    """Contiguous finite derivative-majorant family on the official late plateau.

    The supplied majorants must cover exactly derivative degrees ``0..N`` for
    one spatial window.  Every member is independently passed through
    :func:`section10_endpoint_localization_transfer`, so transition-window or
    non-``t=1`` hypotheses fail closed.

    ``uniform_endpoint_tail_bound`` and ``uniform_cauchy_bound`` take the
    maximum over the finite ladder.  Thus, conditional on the independently
    proved per-degree hypotheses, one obtains one common modulus controlling
    all derivative degrees through ``N`` on the chosen compact spatial window.
    """

    majorants: tuple[EndpointPowerLawMajorant, ...]

    def __init__(self, majorants: Iterable[EndpointPowerLawMajorant]):
        values = tuple(majorants)
        if not values:
            raise ValueError("endpoint majorant ladder must be nonempty")
        if not all(isinstance(item, EndpointPowerLawMajorant) for item in values):
            raise ValueError("every ladder entry must be an EndpointPowerLawMajorant")

        ordered = tuple(sorted(values, key=lambda item: item.derivative_degree))
        degrees = tuple(item.derivative_degree for item in ordered)
        expected = tuple(range(len(ordered)))
        if degrees != expected:
            raise ValueError("endpoint majorant ladder must cover each derivative degree 0..N exactly once")

        windows = {item.spatial_window for item in ordered}
        if len(windows) != 1:
            raise ValueError("endpoint majorant ladder must use one spatial window")

        # This call enforces the exact Section 10 endpoint and the official
        # t>=3/4 plateau for every degree; no sampled interval widening occurs.
        for item in ordered:
            section10_endpoint_localization_transfer(item)

        object.__setattr__(self, "majorants", ordered)

    @property
    def max_degree(self) -> int:
        return len(self.majorants) - 1

    @property
    def spatial_window(self) -> int:
        return self.majorants[0].spatial_window

    @property
    def endpoint(self) -> float:
        return SECTION10_ENDPOINT

    @property
    def common_valid_from(self) -> float:
        return max(item.valid_from for item in self.majorants)

    @property
    def certified(self) -> bool:
        if self.common_valid_from < LATE_START or self.common_valid_from >= self.endpoint:
            return False
        return all(certificate.certified for certificate in self.transfer_certificates())

    def transfer_certificates(self) -> tuple[EndpointLocalizationTransferCertificate, ...]:
        return tuple(section10_endpoint_localization_transfer(item) for item in self.majorants)

    def _common_time(self, value: object, name: str) -> float:
        if isinstance(value, bool):
            raise ValueError(f"{name} must be a finite real number")
        try:
            time = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{name} must be a finite real number") from exc
        if not math.isfinite(time):
            raise ValueError(f"{name} must be a finite real number")
        if time < self.common_valid_from or time >= self.endpoint:
            raise ValueError(
                f"{name} must satisfy common_valid_from <= {name} < endpoint"
            )
        return time

    def endpoint_tail_bounds(self, t: object) -> tuple[float, ...]:
        """Per-degree endpoint tail budgets at one common pre-endpoint time."""
        time = self._common_time(t, "t")
        return tuple(item.endpoint_tail_bound(time) for item in self.majorants)

    def uniform_endpoint_tail_bound(self, t: object) -> float:
        """One finite-order endpoint modulus for all degrees ``0..N``."""
        return max(self.endpoint_tail_bounds(t))

    def cauchy_certificates(
        self, t0: object, t1: object
    ) -> tuple[UniformCauchyCertificate, ...]:
        """Return the per-degree locally-uniform Cauchy implications."""
        start = self._common_time(t0, "t0")
        end = self._common_time(t1, "t1")
        if end < start:
            raise ValueError("Cauchy ladder requires t0 <= t1")
        return tuple(item.cauchy_certificate(start, end) for item in self.majorants)

    def uniform_cauchy_bound(self, t0: object, t1: object) -> float:
        """Maximum interval Cauchy budget across derivative degrees ``0..N``."""
        certificates = self.cauchy_certificates(t0, t1)
        return max(item.interval_upper_bound for item in certificates)
