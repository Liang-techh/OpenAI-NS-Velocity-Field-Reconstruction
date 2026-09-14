"""Pinned scalar tail-order arithmetic for the Section 5 SlowBorel diagonal.

The official ``NavierStokes/DiagonalScale.lean`` theorem
``exists_uniform_tail_order`` turns an increasing gain ``g(j) -> +infinity``
into an arbitrary target power for the *scalar dyadic tail majorant*.  In the
Section 5 specialization used by ``SlowBorelBase.exists_admissibleScales``,

    g(j) = 2 * h * j,

so the theorem's target condition at tail index ``J`` is exactly

    h * (J + 1) - L >= N.

This module makes only that integer-order selection executable, using exact
``Fraction`` arithmetic for the repository's binary64 ``h`` and for ``L,N``.
It does not manufacture the missing all-order hierarchy, the infinite cutoff
schedule, or the PDE residual estimate.  Those hypotheses must exist before
the Lean scalar-tail implication can be applied to the actual background.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .coordinates import validate_h
from .background_truncation_residual import slow_order_exact


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_DIAGONAL_THEOREM = "NavierStokes.DiagonalScale.exists_uniform_tail_order"
PINNED_SLOWBOREL_THEOREM = "NavierStokes.SlowBorelBase.exists_admissibleScales"


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _exact_fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer or Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(f"{name} must be an exact integer or Fraction; floats are forbidden")


def _ceil_fraction(value: Fraction) -> int:
    """Exact ceiling of a Fraction with a positive denominator."""
    return -(-value.numerator // value.denominator)


@dataclass(frozen=True)
class PinnedUniformTailOrderCertificate:
    """Exact integer witness for the pinned scalar target-power inequality.

    ``order`` is the smallest admissible tail index not below
    ``minimum_order`` for which

        g(order + 1) / 2 - derivative_loss >= target_power,

    with the Section 5 gain ``g(j)=2*h*j``.  This is precisely the numerical
    order-selection premise used by ``exists_uniform_tail_order`` after the
    all-order diagonal schedule and scalar tail hypotheses are available.

    The certificate deliberately does *not* assert that those upstream
    hypotheses hold for the reconstructed coefficient hierarchy.
    """

    h: float
    derivative_loss: Fraction | int
    target_power: Fraction | int
    minimum_order: int
    order: int
    formal_revision: str = PINNED_FORMAL_REVISION
    diagonal_theorem: str = PINNED_DIAGONAL_THEOREM
    slowborel_theorem: str = PINNED_SLOWBOREL_THEOREM

    def __post_init__(self) -> None:
        h = validate_h(self.h)
        loss = _exact_fraction(self.derivative_loss, "derivative_loss")
        target = _exact_fraction(self.target_power, "target_power")
        minimum = _nonnegative_int(self.minimum_order, "minimum_order")
        order = _nonnegative_int(self.order, "order")
        if loss < 0:
            raise ValueError("derivative_loss must be nonnegative")
        if target <= 0:
            raise ValueError("target_power must be positive")
        if order < minimum:
            raise ValueError("order must be at least minimum_order")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source revision")
        if self.diagonal_theorem != PINNED_DIAGONAL_THEOREM:
            raise ValueError("diagonal_theorem does not match the pinned theorem")
        if self.slowborel_theorem != PINNED_SLOWBOREL_THEOREM:
            raise ValueError("slowborel_theorem does not match the pinned SlowBorel theorem")

        exact_h = Fraction.from_float(h)
        tail_power = exact_h * (order + 1) - loss
        if tail_power < target:
            raise ValueError("selected order does not reach the requested uniform scalar tail power")

        needed = _ceil_fraction((target + loss) / exact_h) - 1
        expected = max(minimum, max(0, needed))
        if order != expected:
            raise ValueError("order is not the pinned minimal target order above minimum_order")

        # Recheck the same condition in the theorem's displayed gain form.
        gain_next = slow_order_exact(h, order + 1)
        if gain_next != 2 * exact_h * (order + 1):
            raise RuntimeError("slow-order normalization drifted from g(j)=2*h*j")
        if gain_next < 2 * (target + loss):
            raise ValueError("pinned DiagonalScale gain condition is not satisfied")

        object.__setattr__(self, "h", h)
        object.__setattr__(self, "derivative_loss", loss)
        object.__setattr__(self, "target_power", target)
        object.__setattr__(self, "minimum_order", minimum)
        object.__setattr__(self, "order", order)

    @property
    def h_exact(self) -> Fraction:
        return Fraction.from_float(self.h)

    @property
    def gain_next_exact(self) -> Fraction:
        """Return ``g(J+1)=2*h*(J+1)`` exactly for the runtime h."""
        return slow_order_exact(self.h, self.order + 1)

    @property
    def tail_power_exact(self) -> Fraction:
        """Return ``g(J+1)/2-L = h*(J+1)-L`` exactly."""
        return self.gain_next_exact / 2 - self.derivative_loss

    @property
    def target_margin_exact(self) -> Fraction:
        return self.tail_power_exact - self.target_power

    @property
    def dyadic_prefactor_exact(self) -> Fraction:
        """The theorem's exact scalar prefactor ``2^-J``."""
        return Fraction(1, 2) ** self.order

    @property
    def predecessor_fails_target_when_not_minimum_limited(self) -> bool:
        """Check minimality when the requested target, rather than m, fixes J."""
        if self.order == 0 or self.order == self.minimum_order:
            return False
        previous_power = self.h_exact * self.order - self.derivative_loss
        return previous_power < self.target_power

    @property
    def pinned_scalar_tail_arithmetic(self) -> bool:
        return True

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def pde_residual_tail_verified(self) -> bool:
        return False

    @property
    def all_jets_flat(self) -> bool:
        return False

    @property
    def super_algebraic(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False


def certify_pinned_uniform_tail_order(
    h: float,
    derivative_loss: Fraction | int,
    target_power: Fraction | int,
    *,
    minimum_order: int = 0,
) -> PinnedUniformTailOrderCertificate:
    """Select the exact minimal ``J`` for the pinned scalar tail target.

    Floats are forbidden for ``derivative_loss`` and ``target_power`` so a
    rounded near-inequality cannot become a theorem witness.  ``h`` follows the
    repository's existing runtime convention and is interpreted exactly as its
    validated binary64 value, matching ``slow_order_exact``.
    """

    h = validate_h(h)
    loss = _exact_fraction(derivative_loss, "derivative_loss")
    target = _exact_fraction(target_power, "target_power")
    minimum = _nonnegative_int(minimum_order, "minimum_order")
    if loss < 0:
        raise ValueError("derivative_loss must be nonnegative")
    if target <= 0:
        raise ValueError("target_power must be positive")

    exact_h = Fraction.from_float(h)
    needed = _ceil_fraction((target + loss) / exact_h) - 1
    order = max(minimum, max(0, needed))
    return PinnedUniformTailOrderCertificate(
        h=h,
        derivative_loss=loss,
        target_power=target,
        minimum_order=minimum,
        order=order,
    )
