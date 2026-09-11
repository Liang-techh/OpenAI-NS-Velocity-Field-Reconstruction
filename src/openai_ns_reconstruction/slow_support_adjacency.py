"""Section 6 slow-support bridge to the official ``SlotColoring.Adj`` relation.

The landed product cutoff is expressed in normalized chart coordinates
``(R,Z,T)`` on a mesh ``S_*^-3 = ell^-6``.  The pinned official Lean file
``NavierStokes/SlotColoring.lean`` colors *physical* boxes whose axis widths are

    Q^(1/2) ell^-6,  Q^(1/2-h) ell^-6,  Q ell^-6.

These are exactly the physical images of the normalized product mesh under the
landed :class:`~openai_ns_reconstruction.charts.DyadicChart`.  This module makes
that scaling bridge executable and certifies an ``Adj`` witness whenever a
common physical point lies in two landed slow-support enclosures.

No base-field values enter this calculation.  It therefore closes a geometric
formal-structure gap only; it does not instantiate the Proposition 5.5
background, prove the Section 7 uniform base-field estimates, or construct a
paper-exact oscillatory wave.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .charts import DyadicChart
from .coordinates import solve_q_from_tau, validate_h
from .slow_labels import SlowLabel
from .slot_geometry import ExplicitSlotSystem, PairSlotCertificate


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def axis_exponents(h: float) -> tuple[float, float, float]:
    """Physical exponents corresponding to Lean ``axisExponent D``.

    The similarity chart has ``D=1/2-h``, hence the physical scales of
    ``R,Z,T`` are ``Q^(1/2), Q^D, Q`` respectively.
    """
    h = validate_h(h)
    return 0.5, 0.5 - h, 1.0


def physical_mesh_widths(label: SlowLabel, *, h: float = 0.005) -> tuple[float, float, float]:
    """Return the three official physical mesh widths for ``label``.

    This is the executable specialization of ``SlotColoring.spacing`` /
    ``width`` to the manuscript chart exponents.  A zero/overflowed binary64
    width is rejected rather than silently weakening the support certificate.
    """
    h = validate_h(h)
    chart = DyadicChart(label.ell, h)
    mesh = chart.S_star ** -3
    widths = tuple((chart.Q ** exponent) * mesh for exponent in axis_exponents(h))
    if any(not math.isfinite(width) or width <= 0.0 for width in widths):
        raise ArithmeticError("physical mesh width is outside binary64 range")
    return widths


def physical_grid_center(label: SlowLabel, *, h: float = 0.005) -> tuple[float, float, float]:
    """Physical center of the product-grid box identified by ``label.a``."""
    widths = physical_mesh_widths(label, h=h)
    try:
        center = tuple(width * index for width, index in zip(widths, label.a))
    except OverflowError as exc:
        raise ArithmeticError("physical grid center is outside binary64 range") from exc
    if any(not math.isfinite(value) for value in center):
        raise ArithmeticError("physical grid center is outside binary64 range")
    return center


@dataclass(frozen=True)
class PhysicalSlowBox:
    """The official two-mesh physical box used in ``SlotColoring.Adj``.

    ``physicalBox`` in the pinned Lean development has half-width
    ``2 * width`` on every axis.  The actual landed product cutoff has support
    radius one normalized mesh, so it is globally contained in this box and
    leaves one full mesh of enlargement room.
    """

    label: SlowLabel
    h: float = 0.005

    def __post_init__(self) -> None:
        validate_h(self.h)

    @property
    def chart(self) -> DyadicChart:
        return DyadicChart(self.label.ell, self.h)

    @property
    def widths(self) -> tuple[float, float, float]:
        return physical_mesh_widths(self.label, h=self.h)

    @property
    def center(self) -> tuple[float, float, float]:
        return physical_grid_center(self.label, h=self.h)

    def normalized_offsets(self, r: float, z: float, tau: float) -> tuple[float, float, float]:
        """Offsets from the product-grid center measured in mesh units."""
        r, z, tau = (_finite(value, name) for value, name in (
            (r, "r"), (z, "z"), (tau, "tau")
        ))
        if r < 0.0 or tau < 0.0:
            raise ValueError("physical slow point requires r>=0 and tau>=0")
        R, Z, T = self.chart.from_physical_tau(r, z, tau)
        mesh = self.chart.S_star ** -3
        center = tuple(mesh * index for index in self.label.a)
        offsets = tuple((value - c) / mesh for value, c in zip((R, Z, T), center))
        if any(not math.isfinite(value) for value in offsets):
            raise ArithmeticError("normalized support offset is outside binary64 range")
        return offsets

    def contains(self, r: float, z: float, tau: float) -> bool:
        """Check membership in the official two-mesh ``physicalBox``."""
        return all(abs(offset) <= 2.0 for offset in self.normalized_offsets(r, z, tau))

    def contains_landed_cutoff_support_point(self, r: float, z: float, tau: float) -> bool:
        """Check the one-mesh product-cutoff support enclosure at a point."""
        return all(abs(offset) <= 1.0 for offset in self.normalized_offsets(r, z, tau))


@dataclass(frozen=True)
class SlowAdjacencyCertificate:
    """Executable witness of the hypotheses consumed by slot separation."""

    left: SlowLabel
    right: SlowLabel
    q: float
    left_q_ratio: float
    right_q_ratio: float
    product_radius_meshes: float
    left_offsets: tuple[float, float, float]
    right_offsets: tuple[float, float, float]
    slot: PairSlotCertificate

    @property
    def level_gap(self) -> int:
        return abs(self.left.ell - self.right.ell)

    @property
    def certified(self) -> bool:
        return (
            self.left != self.right
            and self.level_gap <= 2
            and 0.5 <= self.left_q_ratio <= 2.0
            and 0.5 <= self.right_q_ratio <= 2.0
            and 0.0 < self.product_radius_meshes <= 2.0
            and all(abs(x) <= self.product_radius_meshes for x in self.left_offsets)
            and all(abs(x) <= self.product_radius_meshes for x in self.right_offsets)
            and self.slot.certified
        )


def certify_slow_support_adjacency(
    left: SlowLabel,
    right: SlowLabel,
    *,
    r: float,
    z: float,
    tau: float,
    h: float = 0.005,
    product_radius_meshes: float = 1.0,
) -> SlowAdjacencyCertificate:
    """Turn a common slow-support point into a ``SlotColoring.Adj`` certificate.

    ``product_radius_meshes=1`` is the closure of the landed
    ``ProductSlowCutoff`` support.  Values up to ``2`` may be used for a
    separately justified fixed enlargement; the official Lean ``physicalBox``
    is exactly the two-mesh envelope.  Larger radii are rejected.

    The common physical point supplies the nonempty intersection required by
    ``Adj``.  The shared physical similarity scale ``q(z,tau)`` is solved once,
    and membership in both dyadic support intervals implies the stronger
    manuscript fact ``|ell-ell'| <= 2`` (hence Lean's required ``<=4``).
    Cross-band product coordinates are checked after converting the *same*
    physical point through each dyadic chart, so different physical mesh
    scalings are not compared as if they were equal.
    """
    if left == right:
        raise ValueError("SlotColoring.Adj requires distinct labels")
    h = validate_h(h)
    radius = _finite(product_radius_meshes, "product_radius_meshes")
    if not 0.0 < radius <= 2.0:
        raise ValueError("product_radius_meshes must lie in (0,2]")

    r, z, tau = (_finite(value, name) for value, name in (
        (r, "r"), (z, "z"), (tau, "tau")
    ))
    if r < 0.0 or tau <= 0.0:
        raise ValueError("interaction witness requires r>=0 and tau>0")

    q = solve_q_from_tau(z, tau, h)
    left_chart = DyadicChart(left.ell, h)
    right_chart = DyadicChart(right.ell, h)
    left_ratio = q / left_chart.Q
    right_ratio = q / right_chart.Q
    if not (0.5 <= left_ratio <= 2.0 and 0.5 <= right_ratio <= 2.0):
        raise ValueError("common point is outside one label's dyadic cutoff support")

    level_gap = abs(left.ell - right.ell)
    # Since Q_ell/Q_m is an exact power of two, two intervals [Q/2,2Q]
    # can share q only when the level gap is <=2.  Check the exact integer
    # consequence rather than accepting a looser sampled/fitted condition.
    if level_gap > 2:
        raise ArithmeticError("dyadic support overlap failed to imply level gap <=2")

    left_box = PhysicalSlowBox(left, h)
    right_box = PhysicalSlowBox(right, h)
    left_offsets = left_box.normalized_offsets(r, z, tau)
    right_offsets = right_box.normalized_offsets(r, z, tau)
    if any(abs(value) > radius for value in left_offsets):
        raise ValueError("common point is outside the left product-support enclosure")
    if any(abs(value) > radius for value in right_offsets):
        raise ValueError("common point is outside the right product-support enclosure")

    # For equal levels, common membership in the official two-mesh boxes gives
    # |a_j-b_j|<=4 exactly, which is the discrete hypothesis used by the mod-5
    # coloring proof.  Check it as an integer invariant before slot_geometry.
    if left.ell == right.ell and any(abs(a - b) > 4 for a, b in zip(left.a, right.a)):
        raise ArithmeticError("same-level physical overlap failed to imply grid gap <=4")

    slot = ExplicitSlotSystem(h).pair_certificate(left, right)
    certificate = SlowAdjacencyCertificate(
        left=left,
        right=right,
        q=q,
        left_q_ratio=left_ratio,
        right_q_ratio=right_ratio,
        product_radius_meshes=radius,
        left_offsets=left_offsets,
        right_offsets=right_offsets,
        slot=slot,
    )
    if not certificate.certified:
        raise ArithmeticError("slow-support adjacency certificate failed")
    return certificate
