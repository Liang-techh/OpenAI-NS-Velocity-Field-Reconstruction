"""Clean-room exact checker for Kokuno NS summation local-band finiteness.

This independently reimplements only the short public local dyadic-band window
argument exposed by the corrected Navier--Stokes reader. It does not import
source-bundle code and it does not prove all-stage convergence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_REPO_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
RESEARCH_STATE_SOURCE_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_COMPONENT = "proof_sources/primary_continuation/summation_realization_body.tex"
SOURCE_COMPONENT_SHA256 = (
    "13d370917df47e46e2d75e0cd3a5d4ab5f2dfd1fd1f31e4eef3fa4413b6a3fa5"
)
SOURCE_CAPTURE_SHA256 = (
    "0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f"
)
BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
DOI = "10.5281/zenodo.22678406"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False


def _q(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, int)):
        raise TypeError(f"{name} must be an exact int or Fraction")
    return Fraction(value)


def dyadic_scale(index: int) -> Fraction:
    """Return Q_index = 2^{-index} exactly."""
    if isinstance(index, bool) or not isinstance(index, int):
        raise TypeError("index must be an int")
    if index >= 0:
        return Fraction(1, 2**index)
    return Fraction(2 ** (-index), 1)


@dataclass(frozen=True)
class LocalBandWindow:
    """Finite conservative dyadic window around one positive local q value."""

    q_local: Fraction
    lower_scale: Fraction
    upper_scale: Fraction
    indices: tuple[int, ...]

    @property
    def index_span(self) -> int:
        if not self.indices:
            return 0
        return self.indices[-1] - self.indices[0]

    @property
    def scales(self) -> tuple[Fraction, ...]:
        return tuple(dyadic_scale(index) for index in self.indices)


def local_band_window(
    q_local: Fraction | int,
    *,
    max_steps: int = 4096,
) -> LocalBandWindow:
    """Enumerate every dyadic Q in the conservative source window [q/4, 4q].

    The public reader first restricts a neighborhood to
    q(point) in (q_local/2, 2 q_local), and an active band obeys
    q(point)/2 <= Q <= 2 q(point). Therefore every active Q is contained in
    the conservative closed interval [q_local/4, 4 q_local]. Since this
    interval has ratio 16, the dyadic index span is at most four.
    """
    q_local = _q(q_local, "q_local")
    if q_local <= 0:
        raise ValueError("q_local must be positive")
    if isinstance(max_steps, bool) or not isinstance(max_steps, int):
        raise TypeError("max_steps must be an int")
    if max_steps <= 0:
        raise ValueError("max_steps must be positive")

    lower = q_local / 4
    upper = 4 * q_local

    index = 0
    scale = Fraction(1)
    steps = 0

    while scale > upper:
        if steps >= max_steps:
            raise RuntimeError("dyadic search exceeded max_steps above window")
        index += 1
        scale /= 2
        steps += 1

    while scale * 2 <= upper:
        if steps >= max_steps:
            raise RuntimeError("dyadic search exceeded max_steps below window")
        index -= 1
        scale *= 2
        steps += 1

    indices: list[int] = []
    while scale >= lower:
        if steps >= max_steps:
            raise RuntimeError("dyadic enumeration exceeded max_steps")
        if scale <= upper:
            indices.append(index)
        index += 1
        scale /= 2
        steps += 1

    if not indices:
        raise AssertionError("a positive ratio-16 window must contain a dyadic scale")
    span = indices[-1] - indices[0]
    if span > 4:
        raise AssertionError("dyadic source window exceeded the source span bound")

    return LocalBandWindow(q_local, lower, upper, tuple(indices))


def active_band_implication(
    q_local: Fraction | int,
    q_point: Fraction | int,
    band_scale: Fraction | int,
) -> bool:
    """Check the exact local hypotheses and their conservative-window conclusion."""
    q_local = _q(q_local, "q_local")
    q_point = _q(q_point, "q_point")
    band_scale = _q(band_scale, "band_scale")
    if q_local <= 0 or q_point <= 0 or band_scale <= 0:
        raise ValueError("q_local, q_point, and band_scale must be positive")

    neighborhood = q_local / 2 < q_point < 2 * q_local
    active_band = q_point / 2 <= band_scale <= 2 * q_point
    if not (neighborhood and active_band):
        return False

    return q_local / 4 <= band_scale <= 4 * q_local


def verify_dyadic_indices(window: LocalBandWindow) -> bool:
    """Verify exact ordering, bounds, completeness by adjacency, and span <= 4."""
    if window.q_local <= 0:
        return False
    if window.lower_scale != window.q_local / 4:
        return False
    if window.upper_scale != 4 * window.q_local:
        return False
    if not window.indices:
        return False
    if tuple(sorted(window.indices)) != window.indices:
        return False
    if len(set(window.indices)) != len(window.indices):
        return False
    scales = window.scales
    if any(not (window.lower_scale <= q <= window.upper_scale) for q in scales):
        return False
    if any(b != a + 1 for a, b in zip(window.indices, window.indices[1:])):
        return False
    first = window.indices[0]
    last = window.indices[-1]
    if window.lower_scale <= dyadic_scale(first - 1) <= window.upper_scale:
        return False
    if window.lower_scale <= dyadic_scale(last + 1) <= window.upper_scale:
        return False
    return window.index_span <= 4


def every_active_band_is_in_window(
    q_local: Fraction | int,
    samples: Iterable[tuple[Fraction | int, Fraction | int]],
) -> bool:
    """Finite regression helper: all supplied valid active samples satisfy the implication."""
    return all(
        active_band_implication(q_local, q_point, band_scale)
        for q_point, band_scale in samples
    )
