"""Typed Section 6 slow-label bridge into the Section 7 phase data.

The paper chooses active slow boxes in Eq. (6.8), duplicates each box with two
sign labels, and freezes the order-zero base jet at a representative before
constructing Eq. (7.2).  This module makes that dependency explicit without
inventing the still-missing squared partition or a surrogate background.

What is executable here:

* the normalized active-shell test in Eq. (6.8), including ``q/Q`` and ``X``;
* the paper mesh scale ``S_*^-3`` and a fail-closed support *enclosure*;
* the label ``gamma=(ell,a,sigma)`` with the two signs sharing one slow box;
* a provider protocol for the actual normalized base jet ``F,G``;
* freezing that jet at the representative and deriving ``g0``, the growth
  frame, ``L_s=2 r0/c_i``, ``k=ceil(epsilon^-1/2)`` and ``B_s``.

The actual product squared partition, Lemma 6.1 rectangle centers/separation,
and the paper-exact Proposition 5.5 background are not constructed here.
Consequently this module is ``formal-structure`` until those upstream objects
supply the provider and support certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Protocol

import numpy as np

from .charts import DyadicChart
from .coordinates import solve_q_from_tau
from .phase import ReferenceFrame, TangentialBaseJet, reference_frame
from .phase_box_certificate import RepresentativePhaseData


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _index3(value, name: str) -> tuple[int, int, int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError(f"{name} must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, int) for x in value):
        raise ValueError(f"{name} must be a length-three integer tuple")
    return value


@dataclass(frozen=True)
class SlowLabel:
    """Paper label ``gamma=(ell,a,sigma)`` from Eq. (6.8)."""

    ell: int
    a: tuple[int, int, int]
    sigma: int

    def __post_init__(self) -> None:
        if isinstance(self.ell, bool) or not isinstance(self.ell, int) or not 1 <= self.ell <= 1000:
            raise ValueError("binary64 labels require integer 1<=ell<=1000")
        _index3(self.a, "a")
        if self.sigma not in (-1, 1):
            raise ValueError("sigma must be +1 or -1")

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        """The sign-free slow box ``beta=(ell,a)`` used by both wave families."""
        return self.ell, self.a

    def opposite(self) -> "SlowLabel":
        return SlowLabel(self.ell, self.a, -self.sigma)


@dataclass(frozen=True)
class SlowBoxEnclosure:
    """Verified mesh-scale enclosure for a labelled slow support.

    Eq. (6.8) uses a product partition whose supports extend by at most one
    mesh length ``S_*^-3`` from a grid point in each slow coordinate.  This
    object verifies only that geometric enclosure.  It deliberately does not
    claim that a particular smooth cutoff has been constructed inside it.
    """

    chart: DyadicChart
    label: SlowLabel
    center: tuple[float, float, float]
    half_width: tuple[float, float, float]

    def __post_init__(self) -> None:
        if self.label.ell != self.chart.ell:
            raise ValueError("label band must equal chart band")
        if len(self.center) != 3 or len(self.half_width) != 3:
            raise ValueError("center and half_width must have three slow coordinates")
        center = tuple(_finite(x, "center") for x in self.center)
        widths = tuple(_finite(x, "half_width") for x in self.half_width)
        if any(x <= 0.0 for x in widths):
            raise ValueError("half_width entries must be positive")
        if any(x > self.mesh_size for x in widths):
            raise ValueError("support enclosure exceeds the paper S_*^-3 mesh scale")
        if center[0] - widths[0] <= 0.0:
            raise ValueError("slow support enclosure must stay away from R=0")

    @property
    def mesh_size(self) -> float:
        return self.chart.S_star ** -3

    def contains(self, R: float, Z: float, T: float) -> bool:
        point = tuple(_finite(x, "slow point") for x in (R, Z, T))
        return all(abs(x - c) <= w for x, c, w in zip(point, self.center, self.half_width))


@dataclass(frozen=True)
class ActiveSlowRepresentative:
    """A representative certified to lie in the active shell of Eq. (6.8).

    ``T0>0`` is used because the repository's executable similarity solver is
    the strict pre-singular branch.  The paper's one-sided smooth extension at
    ``T=0`` remains an analytic boundary statement, not a floating evaluator.
    """

    chart: DyadicChart
    label: SlowLabel
    R0: float
    Z0: float
    T0: float
    q_big: float
    X_a: float
    X_b: float
    enclosure: SlowBoxEnclosure

    def __post_init__(self) -> None:
        if self.label.ell != self.chart.ell or self.enclosure.label.box_key != self.label.box_key:
            raise ValueError("representative, enclosure, and label must use the same slow box")
        R0, Z0, T0 = (_finite(x, name) for x, name in (
            (self.R0, "R0"), (self.Z0, "Z0"), (self.T0, "T0")
        ))
        q_big = _finite(self.q_big, "q_big")
        X_a, X_b = _finite(self.X_a, "X_a"), _finite(self.X_b, "X_b")
        if R0 <= 0.0 or T0 <= 0.0:
            raise ValueError("executable active representative requires R0>0 and T0>0")
        if q_big <= 0.0 or not 0.0 < X_a <= X_b:
            raise ValueError("need q_big>0 and 0<X_a<=X_b")
        if not self.enclosure.contains(R0, Z0, T0):
            raise ValueError("representative must lie in the supplied slow-support enclosure")
        failed = [name for name, ok in self.active_hypotheses().items() if not ok]
        if failed:
            raise ValueError("Eq. (6.8) active-shell hypotheses failed: " + ", ".join(failed))

    @property
    def q_ratio(self) -> float:
        # In chart variables q=Q*s and z=Q^(1/2-h) Z, tau=Q*T, hence
        # s - Z^2 s^(2h) = T.  Reuse the robust scalar solver on (Z,T).
        return solve_q_from_tau(self.Z0, self.T0, self.chart.h)

    @property
    def q(self) -> float:
        return self.chart.Q * self.q_ratio

    @property
    def X0(self) -> float:
        return self.R0 * self.R0 / (2.0 * self.q_ratio)

    @property
    def slow_point(self) -> np.ndarray:
        out = np.array([self.R0, self.Z0, self.T0], dtype=float)
        out.setflags(write=False)
        return out

    def active_hypotheses(self) -> dict[str, bool]:
        ratio = self.q_ratio
        return {
            "q_positive": self.q > 0.0,
            "q_below_q_big": self.q < self.q_big,
            "q_over_Q_band": 0.5 <= ratio <= 2.0,
            "tau_nonnegative": self.T0 >= 0.0,
            "active_annulus": self.X_a <= self.X0 <= self.X_b,
        }


class TangentialBaseJetProvider(Protocol):
    """Future paper-exact base-field hook in normalized chart coordinates."""

    def tangential_jet(
        self, *, chart: DyadicChart, R: float, Z: float, T: float
    ) -> TangentialBaseJet:
        ...


@dataclass(frozen=True)
class FrozenLabelPhaseData:
    """Discrete label data frozen before derivatives, as required in Sections 6-7."""

    source: ActiveSlowRepresentative
    jet0: TangentialBaseJet
    frame: ReferenceFrame
    representative_phase: RepresentativePhaseData
    carrier_k: int
    pulse_length: float
    B_s: float


def freeze_label_phase_data(
    source: ActiveSlowRepresentative,
    provider: TangentialBaseJetProvider,
    *,
    rectangle_radius: float,
    u_star: float,
) -> FrozenLabelPhaseData:
    """Freeze a label's base jet and derive Eqs. (6.11) and (7.2) data.

    ``rectangle_radius`` is the common ``r0`` whose existence is supplied by
    Lemma 6.1; no default is invented.  ``u_star`` is likewise an upstream
    theorem choice.  The provider is called only at the fixed representative.
    """
    r0 = _finite(rectangle_radius, "rectangle_radius")
    u_star = _finite(u_star, "u_star")
    if r0 <= 0.0 or u_star <= 0.0:
        raise ValueError("rectangle_radius and u_star must be positive")

    jet = provider.tangential_jet(
        chart=source.chart, R=source.R0, Z=source.Z0, T=source.T0
    )
    if not isinstance(jet, TangentialBaseJet):
        raise TypeError("provider must return TangentialBaseJet")

    g0 = np.array([source.R0 * jet.F_R, jet.G_R], dtype=float)
    frame = reference_frame(jet.F, g0)
    pulse_length = 2.0 * r0 / source.chart.c_i
    epsilon = source.chart.epsilon
    carrier_k = int(math.ceil(epsilon ** -0.5))
    B_sq = frame.lambda0 / (
        epsilon * float(carrier_k) ** 2 * (1.0 + u_star * u_star) ** 1.5
    )
    if not math.isfinite(B_sq) or B_sq <= 0.0:
        raise ArithmeticError("Eq. (7.2) produced a nonpositive/nonfinite B_s^2")
    B_s = math.sqrt(B_sq)

    phase_rep = RepresentativePhaseData(
        R0=source.R0,
        FR0=jet.F_R,
        GR0=jet.G_R,
        B=B_s,
        sigma=source.label.sigma,
        u=u_star,
        L=pulse_length,
        orientation=1,
    )
    if not np.allclose(phase_rep.K, frame.K, rtol=0.0, atol=32 * np.finfo(float).eps):
        raise ArithmeticError("frozen phase quarter-turn disagrees with the base growth frame")

    return FrozenLabelPhaseData(
        source=source,
        jet0=jet,
        frame=frame,
        representative_phase=phase_rep,
        carrier_k=carrier_k,
        pulse_length=pulse_length,
        B_s=B_s,
    )
