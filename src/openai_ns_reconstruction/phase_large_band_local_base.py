"""Fail-closed large-band LocalBase admission for Sections 6--7.1.

The executable :mod:`slow_labels` path intentionally uses a binary64
``DyadicChart`` and therefore stops at ``ell<=1000``.  The Section 7 estimates,
however, are eventual large-band statements; the exact scale gate in
:mod:`phase_large_band_scale` can work at arbitrarily large integer ``ell``
without ever forming ``Q=2^-ell``.

This module closes one interface gap between those two layers.  It gives a
large-band slow-label address and an admission record for *already certified*
``PhaseEstimates.LocalBaseBounds`` hypotheses.  Error bounds proportional to
``epsilon^2`` are stored as normalized ratios, e.g.

    sup |F-F0| / (M epsilon^2) <= 1,

rather than as binary64 field magnitudes.  Hence neither ``Q`` nor ``epsilon``
is evaluated, and the gate remains meaningful after either quantity would
underflow numerically.

No samples are accepted as theorem evidence.  This module does not construct
Proposition 5.5 base fields, the Section 6 partition, a representative point,
or an oscillatory wave.  Its status is strictly ``formal-structure`` and
``paper_exact_velocity_available`` remains false.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral

from .phase_large_band_scale import LargeBandPhaseScaleCertificate


_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _index3(value: object) -> tuple[int, int, int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError("a must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in value):
        raise ValueError("a must be a length-three integer tuple")
    return tuple(int(x) for x in value)


def _unit_ratio(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite ratio in [0,1]")
    out = float(value)
    if not math.isfinite(out) or not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must be a finite ratio in [0,1]")
    return out


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class AsymptoticSlowLabel:
    """Section 6 label ``gamma=(ell,a,sigma)`` without binary64 chart limits.

    This is an address only.  Unlike :class:`slow_labels.SlowLabel`, it does not
    imply that a floating chart, cutoff, representative, or support has been
    materialized.
    """

    ell: int
    a: tuple[int, int, int]
    sigma: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "ell", _positive_integer(self.ell, "ell"))
        object.__setattr__(self, "a", _index3(self.a))
        if self.sigma not in (-1, 1):
            raise ValueError("sigma must be +1 or -1")

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.ell, self.a

    def opposite(self) -> "AsymptoticSlowLabel":
        return AsymptoticSlowLabel(self.ell, self.a, -self.sigma)


@dataclass(frozen=True)
class LargeBandLocalBaseWitness:
    """Theorem-facing normalized numeric content of ``LocalBaseBounds``.

    Every ``*_over_M`` field represents a certified supremum divided by ``M``.
    Every ``*_over_M_epsilon2`` field represents a certified actual/reference
    error divided by ``M*epsilon^2``.  Values at most one are exactly the
    quantitative hypotheses in the pinned ``PhaseEstimates.LocalBaseBounds``.

    ``diameter_times_S3`` encodes a certified slow-domain diameter ``d`` as
    ``d*S^3``.  Thus ``<=1`` is the large-band-safe form of ``d<=S^-3``.

    The qualitative booleans and ``provenance`` are admission metadata.  They
    are intentionally fail-closed: numerical samples, fits, and empty source
    descriptions are not admissible theorem evidence.
    """

    label: AsymptoticSlowLabel
    M: float
    diameter_times_S3: float
    second_F0_over_M: float
    second_G0_over_M: float
    first_F0_over_M: float
    value_F_error_over_M_epsilon2: float
    first_F_error_over_M_epsilon2: float
    first_G_error_over_M_epsilon2: float
    radial_F_over_M: float
    axial_F_over_M: float
    axial_G_over_M: float
    evidence_kind: str
    provenance: str
    convex_domain_certified: bool
    differentiability_certified: bool
    representative_membership_certified: bool
    local_base_suprema_certified: bool

    def __post_init__(self) -> None:
        if not isinstance(self.label, AsymptoticSlowLabel):
            raise TypeError("label must be an AsymptoticSlowLabel")
        if isinstance(self.M, bool):
            raise ValueError("M must be finite with M>=1")
        M = float(self.M)
        if not math.isfinite(M) or M < 1.0:
            raise ValueError("M must be finite with M>=1")
        object.__setattr__(self, "M", M)

        ratio_names = (
            "diameter_times_S3",
            "second_F0_over_M",
            "second_G0_over_M",
            "first_F0_over_M",
            "value_F_error_over_M_epsilon2",
            "first_F_error_over_M_epsilon2",
            "first_G_error_over_M_epsilon2",
            "radial_F_over_M",
            "axial_F_over_M",
            "axial_G_over_M",
        )
        for name in ratio_names:
            object.__setattr__(self, name, _unit_ratio(getattr(self, name), name))

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; sampled/fitted evidence is rejected"
            )
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("nonempty theorem provenance is required")
        object.__setattr__(self, "provenance", self.provenance.strip())

        for name in (
            "convex_domain_certified",
            "differentiability_certified",
            "representative_membership_certified",
            "local_base_suprema_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def quantitative_hypotheses_admitted(self) -> bool:
        return True

    @property
    def qualitative_hypotheses_admitted(self) -> bool:
        return True


@dataclass(frozen=True)
class LargeBandLocalBaseAdmission:
    """Bind a theorem LocalBase witness to one exact large-band scale.

    This is the large-band counterpart of the numeric ``UniformLocalBaseBridge``
    but deliberately stops before pointwise field evaluation.  It checks band
    identity and ``M`` consistency and preserves ``epsilon^2`` symbolically via
    its exact base-two exponent.
    """

    scale: LargeBandPhaseScaleCertificate
    witness: LargeBandLocalBaseWitness

    def __post_init__(self) -> None:
        if not isinstance(self.scale, LargeBandPhaseScaleCertificate):
            raise TypeError("scale must be a LargeBandPhaseScaleCertificate")
        if not isinstance(self.witness, LargeBandLocalBaseWitness):
            raise TypeError("witness must be a LargeBandLocalBaseWitness")
        if self.witness.label.ell != self.scale.ell:
            raise ValueError("LocalBase label band must equal the large-band phase scale")
        if self.witness.M != self.scale.M:
            raise ValueError("LocalBase M must match the large-band phase scale exactly")
        if not all(self.scale.scale_hypotheses().values()):
            raise ArithmeticError("large-band phase scale unexpectedly lost its certified hypotheses")

    @property
    def ell(self) -> int:
        return self.scale.ell

    @property
    def S_star(self) -> int:
        return self.scale.S_star

    @property
    def max_slow_box_diameter(self) -> Fraction:
        """Exact paper mesh diameter envelope ``S_*^-3=ell^-6``."""
        return Fraction(1, self.S_star**3)

    @property
    def epsilon_squared_log2(self) -> Fraction:
        """Exact ``log2(epsilon^2)=-2 ell h`` without evaluating epsilon."""
        return 2 * self.scale.epsilon_log2

    def admission_checks(self) -> dict[str, bool]:
        return {
            "band_identity": self.witness.label.ell == self.scale.ell,
            "M_identity": self.witness.M == self.scale.M,
            "phase_scale_hypotheses": all(self.scale.scale_hypotheses().values()),
            "diameter_le_S_minus_3": self.witness.diameter_times_S3 <= 1.0,
            "quantitative_local_base_hypotheses": self.witness.quantitative_hypotheses_admitted,
            "qualitative_local_base_hypotheses": self.witness.qualitative_hypotheses_admitted,
            "theorem_evidence_not_sampled": self.witness.evidence_kind in _ACCEPTED_EVIDENCE,
            "provenance_present": bool(self.witness.provenance),
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
