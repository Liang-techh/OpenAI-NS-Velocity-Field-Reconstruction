"""C-infinity experimental cutoffs, not the paper's recursively chosen schedule."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from .coordinates import _finite


def _transition(s: float) -> tuple[float, float]:
    s = _finite(s,"s")
    if s <= 0.5:
        return 1.0, 0.0
    if s >= 1.0:
        return 0.0, 0.0
    v = 2*s-1
    log_ratio = -1/v + 1/(1-v)
    a = math.exp(-abs(log_ratio))
    value = a/(1+a) if log_ratio >= 0 else 1/(1+a)
    derivative = 0.0 if a == 0 else -2*a/(1+a)**2 * (1/v**2 + 1/(1-v)**2)
    return value, derivative


def standard_cutoff(s: float) -> float:
    """Equals one for s<=1/2 and zero for s>=1, flat at both joins."""
    return _transition(s)[0]


def standard_cutoff_derivative(s: float) -> float:
    return _transition(s)[1]


def standard_cutoff_second_derivative(s: float) -> float:
    """Exact second derivative of the explicit transition representative.

    In the transition collar write ``v=2s-1`` and
    ``L=-1/v+1/(1-v)``.  Since ``f=1/(1+exp(L))``, differentiation gives
    ``f''=f(1-f)((1-2f)(L')^2-L'')``.  The flat pieces, including values
    where the stable exponential representation has underflowed all the way
    to an exact endpoint value, have derivative zero.
    """
    s = _finite(s, "s")
    if s <= 0.5 or s >= 1.0:
        return 0.0
    value, _derivative = _transition(s)
    if value == 0.0 or value == 1.0:
        return 0.0
    v = 2.0 * s - 1.0
    l_prime = 2.0 * (1.0 / v**2 + 1.0 / (1.0 - v) ** 2)
    l_second = -8.0 / v**3 + 8.0 / (1.0 - v) ** 3
    second = value * (1.0 - value) * (
        (1.0 - 2.0 * value) * l_prime**2 - l_second
    )
    if not math.isfinite(second):
        raise ArithmeticError("standard cutoff second derivative is not finite")
    return second


def standard_cutoff_third_derivative(s: float) -> float:
    """Exact third derivative of the explicit transition representative.

    With the notation used by :func:`standard_cutoff_second_derivative`,
    ``f'= -f(1-f)L'`` and direct differentiation gives

    ``f''' = f(1-f)[(-1+6f-6f^2)(L')^3 + 3(1-2f)L'L'' - L''']``.

    This derivative belongs only to the repository's executable C-infinity
    representative.  It is not a pointwise realization of Mathlib's
    noncomputable ``ContDiffBump`` in its transition collar.
    """
    s = _finite(s, "s")
    if s <= 0.5 or s >= 1.0:
        return 0.0
    value, _derivative = _transition(s)
    if value == 0.0 or value == 1.0:
        return 0.0
    v = 2.0 * s - 1.0
    l_prime = 2.0 * (1.0 / v**2 + 1.0 / (1.0 - v) ** 2)
    l_second = -8.0 / v**3 + 8.0 / (1.0 - v) ** 3
    l_third = 48.0 / v**4 + 48.0 / (1.0 - v) ** 4
    third = value * (1.0 - value) * (
        (-1.0 + 6.0 * value - 6.0 * value * value) * l_prime**3
        + 3.0 * (1.0 - 2.0 * value) * l_prime * l_second
        - l_third
    )
    if not math.isfinite(third):
        raise ArithmeticError("standard cutoff third derivative is not finite")
    return third


def smooth_cutoff(s: float) -> float:
    """Compatibility name used by the Section-10 geometry module on main."""
    return standard_cutoff(s)


def smooth_step(s: float) -> float:
    """0 on (-inf,0], 1 on [1,inf); C-infinity transition in between."""
    s = _finite(s, "s")
    if s <= 0:
        return 0.0
    if s >= 1:
        return 1.0
    log_ratio = 1/s - 1/(1-s)
    if log_ratio >= 0:
        e = math.exp(-log_ratio)
        return e/(1+e)
    return 1/(1+math.exp(log_ratio))


@dataclass(frozen=True)
class CompactSpatialCutoff:
    """Axisymmetric C-infinity cutoff in r^2+z^2 with an analytic gradient."""
    radius: float = 2.0

    def __post_init__(self) -> None:
        if _finite(self.radius,"radius") <= 0:
            raise ValueError("radius must be positive")

    def __call__(self, x: float, y: float, z: float, t: float) -> float:
        xyz = [_finite(v,n) for v,n in zip((x,y,z),("x","y","z"))]
        return standard_cutoff(sum((v/self.radius)**2 for v in xyz))

    def gradient(self, x: float, y: float, z: float, t: float) -> np.ndarray:
        xyz = np.array([_finite(v,n) for v,n in zip((x,y,z),("x","y","z"))])
        s = float(np.dot(xyz/self.radius,xyz/self.radius))
        return standard_cutoff_derivative(s)*2*xyz/self.radius**2
