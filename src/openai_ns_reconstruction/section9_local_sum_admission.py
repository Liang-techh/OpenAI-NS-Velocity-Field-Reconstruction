"""Fail-closed local-finiteness admission for the Section 9 sum (9.21).

The paper obtains the locally finite sums in Proposition 9.9 from the common
small-q domain of Lemma 9.7 and the shrinking-cutoff mechanism of Lemma 5.4.
For positive correction stages the only arithmetic used here is

    a_{j+1} >= 2 a_j,

with the fixed cutoff support rule ``chi(s)=0`` for ``s>=1``. Therefore, on a
compact q-strip ``c <= q <= c' < q_big``, stage j can be active only if

    a_1 * 2**(j-1) * c < 1.

This module checks that implication exactly with rational arithmetic once an
external theorem/proof witness supplies the common domain and the *infinite*
doubling schedule hypothesis. It does not construct the corrections, infer
an infinite schedule from a finite prefix, or prove Proposition 9.9.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Rational, Real
import math


_RIGOROUS_KINDS = frozenset({
    "paper-derived",
    "lean-derived",
    "certified-numerical",
    "rigorous-external",
})


def _fraction(value: object, name: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive finite rational value")
    if isinstance(value, Fraction):
        out = value
    elif isinstance(value, Integral):
        out = Fraction(int(value), 1)
    elif isinstance(value, Rational):
        out = Fraction(value)
    elif isinstance(value, Real):
        x = float(value)
        if not math.isfinite(x):
            raise ValueError(f"{name} must be finite")
        out = Fraction.from_float(x)
    else:
        raise TypeError(f"{name} must be a real rationalizable value")
    if out <= 0:
        raise ValueError(f"{name} must be positive")
    return out


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _ceil_log2_fraction(value: Fraction) -> int:
    """Return min n>=0 with 2**n >= value, for positive exact ``value``."""
    if value <= 1:
        return 0
    num, den = value.numerator, value.denominator
    n = max(0, num.bit_length() - den.bit_length())
    while (den << n) < num:
        n += 1
    while n > 0 and (den << (n - 1)) >= num:
        n -= 1
    return n


@dataclass(frozen=True)
class Section9LocalSumHypothesis:
    """External theorem input for the Eq. (9.21) local-sum arithmetic.

    Semantically, ``q_big`` is one common domain ``0<q<q_big`` for every
    correction, and ``first_scale`` is ``a_1`` of an *infinite* sequence with
    ``a_{j+1} >= 2 a_j``. ``kind``/``provenance`` identify the independent
    source for those statements; sampled or fitted evidence is intentionally
    not admitted.
    """

    q_big: Fraction | float | int
    first_scale: int
    kind: str
    provenance: str

    def __post_init__(self) -> None:
        q_big = _fraction(self.q_big, "q_big")
        first_scale = _positive_int(self.first_scale, "first_scale")
        if q_big > 1:
            raise ValueError("Lemma 5.4 common domain requires q_big<=1")
        if Fraction(1, first_scale) >= q_big:
            raise ValueError("Lemma 5.4 requires the strict condition 1/a_1 < q_big")
        if self.kind not in _RIGOROUS_KINDS:
            raise ValueError(f"kind must be one of {sorted(_RIGOROUS_KINDS)}")
        provenance = str(self.provenance).strip()
        if not provenance:
            raise ValueError("provenance must be nonempty")
        object.__setattr__(self, "q_big", q_big)
        object.__setattr__(self, "first_scale", first_scale)
        object.__setattr__(self, "provenance", provenance)


@dataclass(frozen=True)
class Section9LocalSumAdmissionCertificate:
    q_big: Fraction
    q_lower: Fraction
    q_upper: Fraction
    first_scale: int
    max_potentially_active_stage: int
    first_guaranteed_inactive_stage: int
    evidence_kind: str
    evidence_provenance: str
    common_domain_arithmetic_verified: bool = True
    local_finiteness_from_doubling_verified: bool = True
    paper_cutoff_support_rule_used: bool = True
    source_theorems_machine_verified: bool = False
    actual_correction_fields_verified: bool = False
    eq_9_21_sum_constructed: bool = False
    proposition_9_9_verified: bool = False
    paper_exact_velocity_available: bool = False


def certify_eq_9_21_local_finiteness(
    hypothesis: Section9LocalSumHypothesis,
    *,
    q_lower: Fraction | float | int,
    q_upper: Fraction | float | int,
) -> Section9LocalSumAdmissionCertificate:
    """Bound the positive stages that can be active on one compact q-strip.

    The returned stage bound follows only from the supplied *infinite* schedule
    theorem input and the paper's fixed support rule. It does not inspect
    correction values and cannot promote a finite list of scales to an
    infinite schedule theorem.
    """
    if not isinstance(hypothesis, Section9LocalSumHypothesis):
        raise TypeError("hypothesis must be a Section9LocalSumHypothesis")
    lower = _fraction(q_lower, "q_lower")
    upper = _fraction(q_upper, "q_upper")
    if not lower < upper:
        raise ValueError("compact Section 9 strip requires q_lower < q_upper")
    if not upper < hypothesis.q_big:
        raise ValueError("compact Section 9 strip requires q_upper < q_big")

    # Stage j can be active only when a_j*q < 1. Since
    # a_j >= a_1*2**(j-1) and q>=q_lower, potential activity implies
    # 2**(j-1) < 1/(a_1*q_lower). The exact number of positive integers j
    # satisfying that strict inequality is ceil(log2(1/(a_1*q_lower))).
    ratio = Fraction(1, hypothesis.first_scale) / lower
    max_active = _ceil_log2_fraction(ratio)
    return Section9LocalSumAdmissionCertificate(
        q_big=hypothesis.q_big,
        q_lower=lower,
        q_upper=upper,
        first_scale=hypothesis.first_scale,
        max_potentially_active_stage=max_active,
        first_guaranteed_inactive_stage=max_active + 1,
        evidence_kind=hypothesis.kind,
        evidence_provenance=hypothesis.provenance,
    )
