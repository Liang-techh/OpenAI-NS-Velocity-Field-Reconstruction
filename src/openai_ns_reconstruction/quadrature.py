"""Cached Gauss-Legendre rules. Numerical integration is not certification."""
from functools import lru_cache
from numbers import Integral
from typing import Callable
import math
import numpy as np


@lru_cache(maxsize=32)
def unit_rule(n: int = 32) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(n, bool) or not isinstance(n, Integral) or not 2 <= n <= 2048:
        raise ValueError("quadrature order must be an integer in [2, 2048]")
    x, w = np.polynomial.legendre.leggauss(int(n))
    x, w = (x + 1) / 2, w / 2
    x.setflags(write=False)
    w.setflags(write=False)
    return x, w


def integrate(fn: Callable[[float], float], a: float, b: float, *, n: int = 32) -> float:
    a, b = float(a), float(b)
    if not math.isfinite(a) or not math.isfinite(b) or b < a:
        raise ValueError("integration endpoints must be finite and ordered")
    x, w = unit_rule(n)
    if a == b:
        return 0.0
    values = np.array([fn(float(a + (b - a) * v)) for v in x], dtype=float)
    if values.shape != x.shape or not np.all(np.isfinite(values)):
        raise ValueError("integrand must return finite scalars")
    return (b - a) * float(w @ values)
