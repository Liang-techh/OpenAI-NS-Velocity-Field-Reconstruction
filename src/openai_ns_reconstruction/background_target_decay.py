"""Finite physical-power targets for certified Section 5 SlowBorel tails.

This module is deliberately downstream of the hierarchy-bound recursive cutoff,
exact retained-recurrence cancellation, and successive residual-extension gates.
It does not manufacture coefficient bounds, recurrence identities, or a new
cutoff schedule.  Instead it makes one missing piece of the finite-prefix
convergence bookkeeping explicit: the first omitted slow order already carries
an exact physical q-exponent

    base_power + 2 * (N + 1) * h.

For a coherent finite chain we can therefore ask, with exact rational arithmetic,
whether each truncation endpoint reaches a requested algebraic q-power.  The
coefficient prefactor is retained in a normalized log-majorant, so a target is
not inferred from an exponent while silently discarding the omitted-tail
coefficient size.

The certificate is intentionally finite and fail-closed.  It does not quantify
over every target power, does not prove a uniform-in-q coefficient envelope,
and therefore does not establish Proposition 5.3, all-jets-flatness, or
super-algebraic convergence.  Those conclusions remain blocked on an actual
all-order hierarchy provider with uniform analytic bounds and exact recurrence
rows.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
import math

from .background_prefix_extension import FiniteCoherentSlowBorelResidualChain
from .background_recurrence_cancellation import CertifiedFullResidualTailMajorant
from .background_truncation_residual import slow_order_exact


def _exact_positive_fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact positive integer or Fraction, not bool")
    if isinstance(value, Fraction):
        result = value
    elif isinstance(value, Integral):
        result = Fraction(int(value), 1)
    else:
        raise TypeError(f"{name} must be an exact integer or Fraction; floats are forbidden")
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _chain_endpoints(
    chain: FiniteCoherentSlowBorelResidualChain,
) -> tuple[CertifiedFullResidualTailMajorant, ...]:
    return (chain.steps[0].previous_residual,) + tuple(
        step.extended_residual for step in chain.steps
    )


@dataclass(frozen=True)
class FinitePhysicalPowerLadderCertificate:
    """Exact target-power admission along one coherent finite truncation chain.

    ``target_exponents[k]`` is the requested physical q-power at endpoint k.
    Targets must be exact integers/Fractions and strictly increase.  Every
    endpoint is rechecked against the physical exponent formula used by the
    first-omitted tail majorant, and the target may not exceed that exponent.

    ``normalized_log_majorants[k]`` bounds the certified residual divided by
    ``q**target_exponents[k]`` at the already-certified q.  It keeps the
    omitted-tail coefficient prefactor visible; it is not a uniform-in-q
    asymptotic constant.
    """

    chain: FiniteCoherentSlowBorelResidualChain
    target_exponents: tuple[Fraction | int, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.chain, FiniteCoherentSlowBorelResidualChain):
            raise TypeError("chain must be a FiniteCoherentSlowBorelResidualChain")
        if not isinstance(self.target_exponents, tuple):
            raise TypeError("target_exponents must be a tuple of exact integers/Fractions")

        endpoints = _chain_endpoints(self.chain)
        if len(self.target_exponents) != len(endpoints):
            raise ValueError(
                "target_exponents must contain exactly one target for each finite chain endpoint"
            )
        targets = tuple(
            _exact_positive_fraction(value, f"target_exponents[{index}]")
            for index, value in enumerate(self.target_exponents)
        )
        if any(right <= left for left, right in zip(targets, targets[1:])):
            raise ValueError("finite physical power targets must be strictly increasing")

        q = endpoints[0].tail.q
        if not 0.0 < q < 1.0:
            raise ValueError("physical power decay requires the certified regime 0 < q < 1")

        previous_exponent: Fraction | None = None
        for index, (residual, target) in enumerate(zip(endpoints, targets)):
            tail = residual.tail
            if tail.q != q:
                raise ValueError("finite physical power ladder changes q")
            expected_exponent = (
                Fraction.from_float(tail.base_power)
                + slow_order_exact(tail.h, residual.first_omitted_order)
            )
            if tail.exponent_exact != expected_exponent:
                raise ValueError(
                    f"endpoint {index} first-omitted exponent is not the physical slow-weight exponent"
                )
            if previous_exponent is not None:
                expected_gain = slow_order_exact(tail.h, 1)
                if tail.exponent_exact - previous_exponent != expected_gain:
                    raise ValueError(
                        "successive physical first-omitted exponents must gain exactly 2*h"
                    )
            if tail.exponent_exact < target:
                raise ValueError(
                    f"endpoint {index} does not attain requested physical q-power {target}"
                )
            if tail.coefficient_l1 < 0.0 or not math.isfinite(tail.coefficient_l1):
                raise ValueError("omitted-tail coefficient prefactor must be finite and nonnegative")
            if math.isnan(tail.log_bound):
                raise ValueError("residual log-majorant must not be NaN")
            previous_exponent = tail.exponent_exact

        object.__setattr__(self, "target_exponents", targets)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def finite_prefix_only(self) -> bool:
        return True

    @property
    def infinite_coherent_family(self) -> bool:
        return False

    @property
    def all_jets_flat(self) -> bool:
        return False

    @property
    def super_algebraic(self) -> bool:
        return False

    @property
    def uniform_in_q(self) -> bool:
        return False

    @property
    def endpoints(self) -> tuple[CertifiedFullResidualTailMajorant, ...]:
        return _chain_endpoints(self.chain)

    @property
    def physical_exponents_exact(self) -> tuple[Fraction, ...]:
        return tuple(endpoint.tail.exponent_exact for endpoint in self.endpoints)

    @property
    def exponent_margins_exact(self) -> tuple[Fraction, ...]:
        return tuple(
            exponent - target
            for exponent, target in zip(self.physical_exponents_exact, self.target_exponents)
        )

    @property
    def normalized_log_majorants(self) -> tuple[float, ...]:
        """Return log bounds for ``residual / q**target`` at the certified q."""
        q = self.endpoints[0].tail.q
        log_q = math.log(q)
        normalized: list[float] = []
        for endpoint, margin in zip(self.endpoints, self.exponent_margins_exact):
            coefficient = endpoint.tail.coefficient_l1
            if coefficient == 0.0:
                normalized.append(-math.inf)
            else:
                normalized.append(math.log(coefficient) + float(margin) * log_q)
        return tuple(normalized)


def certify_finite_physical_power_ladder(
    chain: FiniteCoherentSlowBorelResidualChain,
    target_exponents: tuple[Fraction | int, ...],
) -> FinitePhysicalPowerLadderCertificate:
    """Certify increasing finite q-power targets without an all-order claim."""

    return FinitePhysicalPowerLadderCertificate(
        chain=chain,
        target_exponents=target_exponents,
    )
