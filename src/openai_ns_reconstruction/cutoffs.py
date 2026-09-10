"""C-infinity cutoffs for experiments; not a certified paper cutoff schedule."""
from __future__ import annotations
import math


def smooth_step(s: float) -> float:
    """0 on (-inf,0], 1 on [1,inf); flat to every order at both edges."""
    s = float(s)
    if not math.isfinite(s):
        raise ValueError("cutoff argument must be finite")
    if s <= 0:
        return 0.0
    if s >= 1:
        return 1.0
    log_ratio = 1 / s - 1 / (1 - s)
    if log_ratio >= 0:
        e = math.exp(-log_ratio)
        return e / (1 + e)
    return 1 / (1 + math.exp(log_ratio))


def smooth_cutoff(s: float) -> float:
    """One for s<=1/2, zero for s>=1, with a C-infinity transition."""
    # Avoid cancellation near the zero edge, unlike 1-smooth_step(2*s-1).
    return smooth_step(2 * (1 - float(s)))
