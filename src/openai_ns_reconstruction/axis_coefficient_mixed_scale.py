"""Generic sparse mixed-scale algebra for coefficient eta-jet families.

The wide first-Picard adapters use a few fixed channel widths.  Later Picard
steps need the same arithmetic without deciding those widths in advance.  A
``MixedScaleCoefficient`` stores a sparse immutable map

``(q, p) -> numerator``

for the formal scale ``a(eta)^q * Lambda^-p``.  The numerator is the *full*
eta derivative divided by the common ``a(eta)^q`` factor.  Consequently an eta
derivative is represented by asking the source family for its next jet order;
it is not a derivative of the stored numerator alone.

The amplitude power ``q`` is any nonnegative integer; odd ``q`` channels are
needed for the physical ratio ``F = a phi``.  The algebra never forms ``a^q``
or ``Lambda^-p``.  It preserves every nonzero
channel and performs Decimal arithmetic in a local 96-digit context.  The
family product applies the pinned radial convolution and eta Leibniz rule.
The radial primitives and ``J_r`` are the zero-at-axis operators used by the
pinned coefficient algebra.  Amplitude multiplication must be built as a
family product with a compatible full-derivative Bell family; this module
deliberately provides no bare amplitude-key shift.

This is a formal algebra layer only.  It does not select a SchedulePressure
datum, evaluate a Picard map, certify a norm, or claim convergence to a fixed
point.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import math
from types import MappingProxyType
from typing import Callable, TypeAlias


_DECIMAL_PRECISION = 96

Channel: TypeAlias = tuple[int, int]
ChannelInput: TypeAlias = Mapping[Channel, Decimal] | Iterable[tuple[Channel, Decimal]]
MixedScaleFamily: TypeAlias = Callable[[int, int, float], "MixedScaleCoefficient"]


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _channel_key(value: Channel) -> Channel:
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError("channel keys must be (amplitude_power, inverse_Lambda_power) tuples")
    q, p = value
    if isinstance(q, bool) or not isinstance(q, int) or q < 0:
        raise ValueError("amplitude_power q must be a nonnegative integer")
    if isinstance(p, bool) or not isinstance(p, int) or p < 0:
        raise ValueError("inverse_Lambda_power p must be a nonnegative integer")
    return q, p


def _decimal(value: Decimal | int, name: str) -> Decimal:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a Decimal or integer")
    if isinstance(value, int):
        value = Decimal(value)
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _check_scalar(value: Decimal | int) -> Decimal:
    # Reject float scalars deliberately: this layer must not introduce a
    # binary64 conversion into a theorem-scale channel numerator.
    if isinstance(value, float):
        raise TypeError("mixed-scale scalars must be Decimal or integer, not float")
    return _decimal(value, "scalar")


@dataclass(frozen=True, init=False)
class MixedScaleCoefficient(Mapping[Channel, Decimal]):
    """An immutable sparse map from scale powers to Decimal numerators.

    A key ``(q, p)`` represents the formal channel
    ``a(eta)^q * Lambda^-p``.  Only the numerator is stored.  Zero entries are
    removed, while nonzero entries are never truncated or projected into a
    fixed channel width.
    """

    _items: tuple[tuple[Channel, Decimal], ...] = field(repr=False)

    def __init__(self, channels: ChannelInput | "MixedScaleCoefficient" | None = None) -> None:
        if channels is None:
            entries: Iterable[tuple[Channel, Decimal]] = ()
        elif isinstance(channels, MixedScaleCoefficient):
            entries = channels._items
        elif isinstance(channels, Mapping):
            entries = channels.items()
        else:
            entries = channels

        totals: dict[Channel, Decimal] = {}
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for raw_key, raw_value in entries:
                key = _channel_key(raw_key)
                value = _decimal(raw_value, f"channel {key} numerator")
                if value == 0:
                    continue
                totals[key] = +(totals.get(key, Decimal(0)) + value)
            items = tuple(
                (key, +value)
                for key, value in sorted(totals.items())
                if value != 0
            )
        object.__setattr__(self, "_items", items)

    @classmethod
    def zero(cls) -> "MixedScaleCoefficient":
        """Return the empty coefficient map."""

        return cls()

    @classmethod
    def channel(
        cls,
        q: int,
        p: int,
        numerator: Decimal | int,
    ) -> "MixedScaleCoefficient":
        """Construct one nonzero or zero channel after validating its powers."""

        key = _channel_key((q, p))
        return cls({key: _decimal(numerator, "channel numerator")})

    @property
    def channels(self) -> Mapping[Channel, Decimal]:
        """Return a read-only mapping view of all stored channels."""

        return MappingProxyType(dict(self._items))

    @property
    def support(self) -> tuple[Channel, ...]:
        """Return the sorted nonzero channel keys."""

        return tuple(key for key, _ in self._items)

    def __iter__(self) -> Iterator[Channel]:
        return (key for key, _ in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, key: Channel) -> Decimal:
        for item_key, value in self._items:
            if item_key == key:
                return value
        raise KeyError(key)

    def items(self):
        return self.channels.items()

    def keys(self):
        return self.channels.keys()

    def values(self):
        return self.channels.values()

    def add(self, *others: "MixedScaleCoefficient") -> "MixedScaleCoefficient":
        """Add coefficient maps coefficientwise at Decimal96 precision."""

        totals: dict[Channel, Decimal] = {}
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for coefficient in (self, *others):
                if not isinstance(coefficient, MixedScaleCoefficient):
                    raise TypeError("coefficient addition requires MixedScaleCoefficient operands")
                for key, value in coefficient._items:
                    totals[key] = +(totals.get(key, Decimal(0)) + value)
        return MixedScaleCoefficient(totals)

    def scale(self, scalar: Decimal | int) -> "MixedScaleCoefficient":
        """Multiply every numerator by a finite Decimal or integer scalar."""

        scalar_decimal = _check_scalar(scalar)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return MixedScaleCoefficient(
                {
                    key: +(scalar_decimal * value)
                    for key, value in self._items
                }
            )

    def shift_lambda(self, delta: int) -> "MixedScaleCoefficient":
        """Shift every inverse-Lambda power by ``delta``.

        Lambda is eta-independent, so this ordinary scale shift is valid.  A
        shift that would create a negative inverse power is rejected.  No
        amplitude-power shift is provided because it would be wrong for full
        eta derivatives unless the Bell derivative family were multiplied too.
        """

        if isinstance(delta, bool) or not isinstance(delta, int):
            raise TypeError("Lambda power shift must be an integer")
        shifted: dict[Channel, Decimal] = {}
        for (q, p), value in self._items:
            shifted_power = p + delta
            if shifted_power < 0:
                raise ValueError("Lambda power shift would create a negative inverse power")
            shifted[(q, shifted_power)] = value
        return MixedScaleCoefficient(shifted)

    def multiply(self, other: "MixedScaleCoefficient") -> "MixedScaleCoefficient":
        """Multiply channels by adding both scale powers, without truncation."""

        if not isinstance(other, MixedScaleCoefficient):
            raise TypeError("channel multiplication requires MixedScaleCoefficient")
        totals: dict[Channel, Decimal] = {}
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for (left_q, left_p), left_value in self._items:
                for (right_q, right_p), right_value in other._items:
                    key = (left_q + right_q, left_p + right_p)
                    totals[key] = +(
                        totals.get(key, Decimal(0))
                        + left_value * right_value
                    )
        return MixedScaleCoefficient(totals)

    def __add__(self, other: "MixedScaleCoefficient") -> "MixedScaleCoefficient":
        if not isinstance(other, MixedScaleCoefficient):
            return NotImplemented
        return self.add(other)

    def __radd__(self, other: object) -> "MixedScaleCoefficient":
        if other == 0:
            return self
        if not isinstance(other, MixedScaleCoefficient):
            return NotImplemented
        return other.add(self)

    def __sub__(self, other: "MixedScaleCoefficient") -> "MixedScaleCoefficient":
        if not isinstance(other, MixedScaleCoefficient):
            return NotImplemented
        return self.add(other.scale(-1))

    def __neg__(self) -> "MixedScaleCoefficient":
        return self.scale(-1)

    def __mul__(self, scalar: Decimal | int) -> "MixedScaleCoefficient":
        if isinstance(scalar, (Decimal, int)) and not isinstance(scalar, bool):
            return self.scale(scalar)
        return NotImplemented

    def __rmul__(self, scalar: Decimal | int) -> "MixedScaleCoefficient":
        return self.__mul__(scalar)


def _family_value(
    provider: MixedScaleFamily,
    n: int,
    m: int,
    eta: float,
) -> MixedScaleCoefficient:
    if not callable(provider):
        raise TypeError("mixed-scale family provider must be callable")
    value = provider(n, m, eta)
    if isinstance(value, MixedScaleCoefficient):
        return value
    if isinstance(value, Mapping):
        return MixedScaleCoefficient(value)
    raise TypeError("mixed-scale family provider must return a channel mapping")


def mixed_scale_product(
    left: MixedScaleFamily,
    right: MixedScaleFamily,
) -> MixedScaleFamily:
    """Build the pinned radial/eta product of two mixed-scale families."""

    if not callable(left) or not callable(right):
        raise TypeError("mixed-scale product operands must be callable families")

    def product(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        m = _index(m, "m")
        totals: dict[Channel, Decimal] = {}
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k))
                    term = _family_value(left, i, k, eta).multiply(
                        _family_value(right, j, l, eta)
                    )
                    for key, value in term._items:
                        totals[key] = +(totals.get(key, Decimal(0)) + weight * value)
        return MixedScaleCoefficient(totals)

    return product


def mixed_scale_eta_derivative(source: MixedScaleFamily) -> MixedScaleFamily:
    """Build the full eta derivative family by requesting jet order ``m+1``.

    The source channel numerators already represent full derivatives divided
    by their common amplitude power.  This operation therefore performs the
    exact family index shift and does not differentiate a numerator in place.
    """

    if not callable(source):
        raise TypeError("eta-derivative operand must be a callable family")

    def derivative(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        m = _index(m, "m")
        return _family_value(source, n, m + 1, eta)

    return derivative


def mixed_scale_euler(source: MixedScaleFamily) -> MixedScaleFamily:
    """Build the radial Euler derivative family ``n * source[n,m]``."""

    if not callable(source):
        raise TypeError("Euler-derivative operand must be a callable family")

    def euler(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        m = _index(m, "m")
        return _family_value(source, n, m, eta).scale(n)

    return euler


def mixed_scale_average(source: MixedScaleFamily) -> MixedScaleFamily:
    """Build the pinned radial average family ``source[n,m] / (n+1)``."""

    if not callable(source):
        raise TypeError("average operand must be a callable family")

    def average(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        m = _index(m, "m")
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            scale = +(Decimal(1) / Decimal(n + 1))
        return _family_value(source, n, m, eta).scale(scale)

    return average


def mixed_scale_primitive(source: MixedScaleFamily) -> MixedScaleFamily:
    """Build the zero-at-axis radial primitive ``source[n-1,m] / n``."""

    if not callable(source):
        raise TypeError("primitive operand must be a callable family")

    def primitive(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        m = _index(m, "m")
        if n == 0:
            # The zero-at-axis primitive is structural.  Do not probe the
            # source row here: a triangular provider may depend on this same
            # output row and must remain causally well-founded.
            return MixedScaleCoefficient.zero()
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            scale = +(Decimal(1) / Decimal(n))
        return _family_value(source, n - 1, m, eta).scale(scale)

    return primitive


def _j_radius(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value not in (1, 2):
        raise ValueError("J_r requires r equal to 1 or 2")
    return value


def mixed_scale_j_r(source: MixedScaleFamily, r: int) -> MixedScaleFamily:
    """Build the pinned zero-at-axis ``J_r`` family for ``r = 1`` or ``2``."""

    if not callable(source):
        raise TypeError("J_r operand must be a callable family")
    r = _j_radius(r)

    def inverse(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        m = _index(m, "m")
        if n == 0:
            # The zero-at-axis J_r row is structural.  Avoid probing the
            # source row so triangular providers cannot recurse at n = 0.
            return MixedScaleCoefficient.zero()
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n * (n + r - 1))
            scale = +(Decimal(1) / divisor)
        return _family_value(source, n - 1, m, eta).scale(scale)

    return inverse


__all__ = [
    "Channel",
    "ChannelInput",
    "MixedScaleCoefficient",
    "MixedScaleFamily",
    "mixed_scale_average",
    "mixed_scale_eta_derivative",
    "mixed_scale_euler",
    "mixed_scale_j_r",
    "mixed_scale_primitive",
    "mixed_scale_product",
]
