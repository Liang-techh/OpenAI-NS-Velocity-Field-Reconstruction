"""Fail-closed finite-power admission for the Section 9 flat remainder.

Lemma 9.8 / Proposition 9.9 require the error term ``E_{j,m}`` in Eq. (9.18)
to be flat in the dyadic scale: for every power ``N`` there is a theorem-level
constant ``C_{j,m,N}`` such that, on one common small-``q`` domain,

    E_{j,m}(q) <= C_{j,m,N} q**N.

A pointwise residual check at one supplied ``q`` cannot establish this.  This
module therefore does not evaluate residual samples and does not fit power laws.
Instead it defines an audit/admission interface for externally justified
power-law witnesses and verifies that a *finite* family covers the exact
contiguous powers ``0..N`` on one common stage, derivative order and ``q``
domain.

Passing this gate proves only that the supplied finite witness family is
structurally coherent.  It is not the all-orders statement, does not construct
the genuine Section 9 correction sequence, and does not promote the repository
to paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import operator
from typing import Iterable

from .section9_stage_certificate import CertifiedBoundDatum


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(out)


def _q_upper(value: object) -> float:
    if isinstance(value, bool):
        raise ValueError("q_upper must satisfy 0 < q_upper < 1")
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("q_upper must satisfy 0 < q_upper < 1") from exc
    if not math.isfinite(out) or not 0.0 < out < 1.0:
        raise ValueError("q_upper must satisfy 0 < q_upper < 1")
    return out


@dataclass(frozen=True)
class Section9FlatRemainderPowerWitness:
    """External theorem-level witness for one power of ``E_{j,m}``.

    Semantics: ``constant_bound`` certifies that on the full common domain
    ``0 < q <= q_upper`` the *same* Section 9 flat remainder satisfies

    ``E_{j,m}(q) <= constant_bound.upper_bound * q**power``.

    ``constant_bound`` is reused from :mod:`section9_stage_certificate` so an
    ordinary sampled/fitted estimate cannot enter this interface without an
    allowed rigorous evidence classification and nonempty provenance string.
    This class validates the metadata but cannot prove the external theorem
    named by that provenance.
    """

    stage: int
    derivative_order: int
    power: int
    q_upper: float
    constant_bound: CertifiedBoundDatum
    status: str = "formal-structure"
    uniform_power_law_source_verified: bool = False
    actual_section9_sequence_verified: bool = False
    paper_exact_velocity_available: bool = False

    def __post_init__(self) -> None:
        stage = _natural(self.stage, "stage")
        derivative_order = _natural(self.derivative_order, "derivative_order")
        power = _natural(self.power, "power")
        q_upper = _q_upper(self.q_upper)
        if not isinstance(self.constant_bound, CertifiedBoundDatum):
            raise TypeError("constant_bound must be CertifiedBoundDatum")
        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "derivative_order", derivative_order)
        object.__setattr__(self, "power", power)
        object.__setattr__(self, "q_upper", q_upper)


@dataclass(frozen=True)
class Section9FlatRemainderLadderCertificate:
    """Coherent finite prefix of the all-powers flatness requirement."""

    stage: int
    derivative_order: int
    q_upper: float
    max_power: int
    witnesses: tuple[Section9FlatRemainderPowerWitness, ...]
    status: str = "formal-structure"
    contiguous_finite_power_ladder_verified: bool = True
    common_q_domain_verified: bool = True
    source_theorems_machine_verified: bool = False
    flat_remainder_all_orders_verified: bool = False
    actual_section9_sequence_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def powers(self) -> tuple[int, ...]:
        return tuple(w.power for w in self.witnesses)

    @property
    def formal_ladder_ready(self) -> bool:
        return (
            self.powers == tuple(range(self.max_power + 1))
            and all(w.stage == self.stage for w in self.witnesses)
            and all(w.derivative_order == self.derivative_order for w in self.witnesses)
            and all(w.q_upper == self.q_upper for w in self.witnesses)
            and all(bool(w.constant_bound.provenance.strip()) for w in self.witnesses)
            and not self.flat_remainder_all_orders_verified
            and not self.actual_section9_sequence_verified
            and not self.paper_exact_velocity_available
        )


def certify_flat_remainder_power_ladder(
    witnesses: Iterable[Section9FlatRemainderPowerWitness],
    *,
    max_power: int,
) -> Section9FlatRemainderLadderCertificate:
    """Validate one finite ``N=0..max_power`` flat-remainder witness ladder.

    The powers must be present exactly once and must refer to the same Section 9
    stage, derivative order and small-``q`` domain.  The function intentionally
    does not extrapolate from this finite family to arbitrary ``N``.
    """

    max_power = _natural(max_power, "max_power")
    rows = tuple(witnesses)
    if not rows:
        raise ValueError("flat-remainder ladder must contain at least one witness")
    if not all(isinstance(row, Section9FlatRemainderPowerWitness) for row in rows):
        raise TypeError("every ladder entry must be Section9FlatRemainderPowerWitness")

    by_power: dict[int, Section9FlatRemainderPowerWitness] = {}
    for row in rows:
        if row.power in by_power:
            raise ValueError(f"duplicate flat-remainder power N={row.power}")
        by_power[row.power] = row

    required = set(range(max_power + 1))
    supplied = set(by_power)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            "flat-remainder ladder must cover exactly powers 0..max_power; "
            f"missing={missing}, extra={extra}"
        )

    ordered = tuple(by_power[n] for n in range(max_power + 1))
    first = ordered[0]
    for row in ordered[1:]:
        if row.stage != first.stage:
            raise ValueError("all flat-remainder witnesses must use one Section 9 stage")
        if row.derivative_order != first.derivative_order:
            raise ValueError("all flat-remainder witnesses must use one derivative order")
        if row.q_upper != first.q_upper:
            raise ValueError("all flat-remainder witnesses must use one common q domain")

    return Section9FlatRemainderLadderCertificate(
        stage=first.stage,
        derivative_order=first.derivative_order,
        q_upper=first.q_upper,
        max_power=max_power,
        witnesses=ordered,
    )
