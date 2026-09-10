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
    # log_ratio is already huge before the algebraic derivative can overflow.
    derivative = 0.0 if a == 0 else -2*a/(1+a)**2 * (1/v**2 + 1/(1-v)**2)
    return value, derivative


def standard_cutoff(s: float) -> float:
    """Equals one for s<=1/2 and zero for s>=1, flat at both joins."""
    return _transition(s)[0]


def standard_cutoff_derivative(s: float) -> float:
    return _transition(s)[1]


@dataclass(frozen=True)
class CompactSpatialCutoff:
    """Axisymmetric C-infinity cutoff in r^2+z^2 with an analytic gradient.

    This is an experimental spatial cutoff; it is not a temporal force extension.
    """
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
