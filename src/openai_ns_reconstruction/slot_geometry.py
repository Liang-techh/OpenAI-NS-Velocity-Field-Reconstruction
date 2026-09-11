"""Constructive auxiliary-slot geometry for Section 6 / Lemma 6.1.

This module supplies a reproducible, band-independent family of auxiliary
centres and a *conservative* common rectangle radius.  It follows the explicit
finite colouring and rational-centre construction in the pinned official Lean
files ``SlotColoring.lean`` and ``SlotGeometry.lean``.

The important boundary is provenance: the manuscript only needs the existence
of a finite proper colouring and a common positive radius.  The official Lean
uses ``Fintype.equivFin`` to enumerate its 2250 palette and obtains the radius
by an open-neighbourhood existence argument.  Here we choose an explicit
mixed-radix palette enumeration and an explicit smaller rational radius.  Both
are admissible witnesses, but neither is claimed to be definitionally equal to
the noncomputable/opaque choices in the paper or Lean development.  This is
therefore ``formal-structure`` until the landed slow supports are connected to
the full physical ``Adj`` relation and the paper-exact background is available.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from .charts import DyadicChart
from .coordinates import validate_h
from .slow_labels import SlowLabel


PALETTE_CARD = 2250
COVER_NORM_BOUND = 6
_ORIENTED_RADIUS_DENOMINATOR = 3  # 1 + ||v_r||_inf + ||v_t||_inf.


def native_gap(h: float) -> int:
    """Official ``SlotColoring.nativeGap`` for the manuscript covering.

    For labels whose dyadic levels differ by at most four, the pinned Lean
    theorem ``nativeIndex_gap`` bounds the difference of their covering
    indices by this integer.
    """
    h = validate_h(h)
    growth = 4.0 + math.sqrt(2.0)
    budget = (4.0 * (1.0 + h) * math.log(2.0) + 2.0 * math.log(5.0)) / math.log(growth)
    return int(math.ceil(budget)) + 1


def palette_coordinates(label: SlowLabel) -> tuple[int, tuple[int, int, int], int]:
    """The explicit ``Fin 9 x (Fin 5)^3 x Bool`` colour data.

    Python's remainder by a positive integer is the same nonnegative residue
    convention used by Lean's ``Nat.mod``/``Int.emod`` here.
    """
    residues = tuple(a % 5 for a in label.a)
    sign_bit = 1 if label.sigma > 0 else 0
    return label.ell % 9, residues, sign_bit


def palette_index(label: SlowLabel) -> int:
    """Deterministically enumerate the 2250 palette entries.

    This mixed-radix enumeration is an executable admissible choice.  It is
    not asserted to coincide with Lean's opaque ``Fintype.equivFin`` ordering.
    """
    level, grid, sign = palette_coordinates(label)
    out = level
    for residue in grid:
        out = 5 * out + residue
    out = 2 * out + sign
    if not 0 <= out < PALETTE_CARD:  # defensive: all arithmetic above is exact.
        raise ArithmeticError("palette enumeration escaped Fin 2250")
    return out


def certify_interaction_colors(left: SlowLabel, right: SlowLabel) -> tuple[int, int]:
    """Certify the proper-colouring step for a paper interaction pair.

    The Section-6 interaction graph uses distinct labels with dyadic-level gap
    at most four.  At equal level, overlap of the enlarged product boxes gives
    coordinate grid gaps at most four.  Those are the exact finite hypotheses
    needed by the ``mod 9 / mod 5 / sign`` colouring argument.

    This function deliberately does *not* infer physical slow-box overlap from
    samples.  It only verifies the discrete consequences supplied by that
    analytic/support step.
    """
    if left == right:
        raise ValueError("interaction colours require distinct labels")
    if abs(left.ell - right.ell) > 4:
        raise ValueError("paper interaction certification requires dyadic-level gap <= 4")
    if left.ell == right.ell and any(abs(a - b) > 4 for a, b in zip(left.a, right.a)):
        raise ValueError("same-level interaction certification requires grid gaps <= 4")
    colors = palette_index(left), palette_index(right)
    if colors[0] == colors[1]:
        raise ArithmeticError("proper-colouring invariant failed")
    return colors


def _matmul2(a: tuple[tuple[int, int], tuple[int, int]],
             b: tuple[tuple[int, int], tuple[int, int]]) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0],
         a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0],
         a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )


def covering_power(n: int) -> tuple[tuple[int, int], tuple[int, int]]:
    """Exact integer power of ``J=[[3,1],[1,5]]``."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("covering power must be a nonnegative integer")
    result = ((1, 0), (0, 1))
    base = ((3, 1), (1, 5))
    k = n
    while k:
        if k & 1:
            result = _matmul2(result, base)
        base = _matmul2(base, base)
        k >>= 1
    return result


def covering_inf_norm(n: int) -> int:
    """Exact induced infinity norm of ``J^n``."""
    matrix = covering_power(n)
    return max(sum(abs(v) for v in row) for row in matrix)


def _distance_to_integer(x: Fraction) -> Fraction:
    floor = x.numerator // x.denominator
    frac = x - floor
    return min(frac, 1 - frac)


def torus_sup_distance(x: tuple[Fraction, Fraction], y: tuple[Fraction, Fraction]) -> Fraction:
    """Exact quotient distance ``min_{k in Z^2} ||x-y-k||_inf``."""
    return max(_distance_to_integer(a - b) for a, b in zip(x, y))


@dataclass(frozen=True)
class PairSlotCertificate:
    """Exact finite certificate for one pair of interacting labels."""

    lower_label: SlowLabel
    upper_label: SlowLabel
    lower_color: int
    upper_color: int
    covering_gap: int
    center_torus_distance: Fraction
    universal_center_gap: Fraction
    separation_perturbation_bound: Fraction
    injection_perturbation_bound: Fraction
    rectangle_radius: Fraction

    @property
    def certified(self) -> bool:
        return (
            self.lower_color != self.upper_color
            and 0 <= self.covering_gap
            and self.center_torus_distance >= self.universal_center_gap > 0
            and self.separation_perturbation_bound < self.universal_center_gap
            and self.injection_perturbation_bound < 1
            and self.rectangle_radius > 0
        )


@dataclass(frozen=True)
class ExplicitSlotSystem:
    """Band-independent Section-6 centres and a common padded radius witness."""

    h: float = 0.005

    def __post_init__(self) -> None:
        validate_h(self.h)

    @property
    def gap(self) -> int:
        return native_gap(self.h)

    @property
    def denominator(self) -> int:
        # ``SlotGeometry.denominator m D = (m+1) * 6^D``.
        return (PALETTE_CARD + 1) * COVER_NORM_BOUND ** self.gap

    @property
    def universal_center_gap(self) -> Fraction:
        """A uniform lower bound for every forbidden rational-centre coincidence."""
        return Fraction(1, self.denominator)

    @property
    def standard_radius(self) -> Fraction:
        """Conservative radius for the sup-norm rectangles in ``SlotGeometry``.

        The two terms respectively leave a factor-two margin for quotient
        injectivity and for separation after any covering gap ``<= D``.
        """
        six_d = COVER_NORM_BOUND ** self.gap
        injective = Fraction(1, 8 * six_d)
        separated = Fraction(1, 4 * self.denominator * (six_d + 1))
        return min(injective, separated)

    @property
    def rectangle_radius(self) -> Fraction:
        """Common coordinate ``r0`` for manuscript-oriented rectangles.

        The manuscript axes are ``v_r=(1,-beta)``, ``v_t=(beta,1)`` with
        ``beta=sqrt(2)-1``, hence both have sup norm one.  The official Lean
        inclusion ``orientedRectangle_subset`` therefore divides the standard
        radius by ``1+1+1=3``.
        """
        return self.standard_radius / _ORIENTED_RADIUS_DENOMINATOR

    @property
    def rectangle_radius_float(self) -> float:
        value = float(self.rectangle_radius)
        if not math.isfinite(value) or value <= 0.0:
            raise ArithmeticError("common rectangle radius is outside binary64 range")
        return value

    def center_for_color(self, color: int) -> tuple[Fraction, Fraction]:
        if isinstance(color, bool) or not isinstance(color, int) or not 0 <= color < PALETTE_CARD:
            raise ValueError("color must be an integer in [0,2250)")
        return Fraction(color + 1, self.denominator), Fraction(0, 1)

    def center(self, label: SlowLabel) -> tuple[Fraction, Fraction]:
        return self.center_for_color(palette_index(label))

    def pair_certificate(self, left: SlowLabel, right: SlowLabel) -> PairSlotCertificate:
        """Certify common-slot separation for one paper-interaction pair.

        The labels are first checked against the discrete consequences of the
        enlarged-support interaction relation.  Their actual covering indices
        are then obtained from the landed ``DyadicChart`` implementation.
        No base-field or phase samples enter this certificate.
        """
        left_color, right_color = certify_interaction_colors(left, right)
        left_level = DyadicChart(left.ell, self.h).covering_index
        right_level = DyadicChart(right.ell, self.h).covering_index
        if left_level <= right_level:
            lo, hi = left, right
            lo_color, hi_color = left_color, right_color
            lo_level, hi_level = left_level, right_level
        else:
            lo, hi = right, left
            lo_color, hi_color = right_color, left_color
            lo_level, hi_level = right_level, left_level
        gap = hi_level - lo_level
        if gap > self.gap:
            raise ArithmeticError("landed chart covering indices violate nativeGap certificate")

        matrix = covering_power(gap)
        c_lo = self.center_for_color(lo_color)
        c_hi = self.center_for_color(hi_color)
        covered = (
            matrix[0][0] * c_lo[0] + matrix[0][1] * c_lo[1],
            matrix[1][0] * c_lo[0] + matrix[1][1] * c_lo[1],
        )
        center_distance = torus_sup_distance(covered, c_hi)
        if center_distance < self.universal_center_gap:
            raise ArithmeticError("explicit rational centers lost their universal lattice gap")

        norm = covering_inf_norm(gap)
        standard = self.standard_radius
        separation_error = 2 * standard * (norm + 1)
        injection_error = 4 * standard * norm
        cert = PairSlotCertificate(
            lower_label=lo,
            upper_label=hi,
            lower_color=lo_color,
            upper_color=hi_color,
            covering_gap=gap,
            center_torus_distance=center_distance,
            universal_center_gap=self.universal_center_gap,
            separation_perturbation_bound=separation_error,
            injection_perturbation_bound=injection_error,
            rectangle_radius=self.rectangle_radius,
        )
        if not cert.certified:
            raise ArithmeticError("common slot separation certificate failed")
        return cert
