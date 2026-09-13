"""Exact retained-recurrence certificates for the Section 5 slow expansion.

``background_truncation_residual`` already exposes the exact finite recurrence /
truncation split and an omitted-tail majorant, but deliberately refuses to
promote a numerically small retained recurrence to zero.  This module supplies
one stricter bridge: retained recurrence cancellation must be represented by an
*exact formal identity* with hierarchy provenance before the omitted-tail
majorant may be interpreted as a bound for the full finite residual.

The formal algebra uses :class:`fractions.Fraction` coefficients only.  Floats
are rejected, so there is no tolerance, sampled residual, or ``isclose`` gate in
the cancellation decision.  The current layer remains finite-prefix
``formal-structure`` until actual hierarchy constructors populate these exact
identities and the all-order argument is closed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Iterable, Sequence

from .background_truncation_residual import SlowTailMajorant, omitted_slow_tail_majorant


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


def _exact_fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer or Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(f"{name} must be an exact integer or Fraction; floats are forbidden")


@dataclass(frozen=True)
class FormalAtomTerm:
    """One exact coefficient multiplying an opaque formal hierarchy atom."""

    atom: str
    coefficient: Fraction | int

    def __post_init__(self) -> None:
        object.__setattr__(self, "atom", _nonempty_text(self.atom, "atom"))
        object.__setattr__(
            self,
            "coefficient",
            _exact_fraction(self.coefficient, "coefficient"),
        )


def _canonical_expression(terms: Iterable[FormalAtomTerm]) -> tuple[FormalAtomTerm, ...]:
    totals: dict[str, Fraction] = {}
    for term in terms:
        if not isinstance(term, FormalAtomTerm):
            raise TypeError("formal expressions must contain FormalAtomTerm values")
        totals[term.atom] = totals.get(term.atom, Fraction(0, 1)) + term.coefficient
    return tuple(
        FormalAtomTerm(atom, coefficient)
        for atom, coefficient in sorted(totals.items())
        if coefficient != 0
    )


def _negated(terms: Iterable[FormalAtomTerm]) -> tuple[FormalAtomTerm, ...]:
    return tuple(FormalAtomTerm(term.atom, -term.coefficient) for term in terms)


@dataclass(frozen=True)
class RetainedRecurrenceDependency:
    """Hierarchy coefficient artifact used by one exact retained identity."""

    coefficient_order: int
    artifact: str
    provider: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "coefficient_order",
            _nonnegative_int(self.coefficient_order, "coefficient_order"),
        )
        object.__setattr__(self, "artifact", _nonempty_text(self.artifact, "artifact"))
        object.__setattr__(self, "provider", _nonempty_text(self.provider, "provider"))


@dataclass(frozen=True)
class ExactRetainedRecurrenceIdentity:
    """Exact function-level form of ``L_n + K_n - previous(A)_n = 0``.

    ``linear_terms``, ``pair_terms`` and ``previous_shifted_terms`` are formal
    atoms with exact rational coefficients.  Construction succeeds only if
    their canonical residual is identically empty.  Consequently a value such
    as ``1/10**30`` is *not* zero and cannot pass this gate.
    """

    order: int
    hierarchy_id: str
    source_revision: str
    coefficient_state_id: str
    identity_source: str
    theorem_name: str
    dependencies: tuple[RetainedRecurrenceDependency, ...]
    linear_terms: tuple[FormalAtomTerm, ...]
    pair_terms: tuple[FormalAtomTerm, ...]
    previous_shifted_terms: tuple[FormalAtomTerm, ...]

    def __post_init__(self) -> None:
        order = _nonnegative_int(self.order, "order")
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        coefficient_state_id = _nonempty_text(self.coefficient_state_id, "coefficient_state_id")
        identity_source = _nonempty_text(self.identity_source, "identity_source")
        theorem_name = _nonempty_text(self.theorem_name, "theorem_name")
        dependencies = tuple(self.dependencies)
        if not dependencies:
            raise ValueError("retained recurrence identity must carry hierarchy dependencies")
        if not all(isinstance(dep, RetainedRecurrenceDependency) for dep in dependencies):
            raise TypeError("dependencies must contain RetainedRecurrenceDependency values")
        dep_keys = [(dep.coefficient_order, dep.artifact, dep.provider) for dep in dependencies]
        if len(dep_keys) != len(set(dep_keys)):
            raise ValueError("retained recurrence dependencies must be unique")
        if any(dep.coefficient_order > order for dep in dependencies):
            raise ValueError("retained recurrence cannot depend on a future coefficient order")
        if not any(dep.coefficient_order == order for dep in dependencies):
            raise ValueError("retained recurrence must bind its own coefficient order")

        linear = _canonical_expression(tuple(self.linear_terms))
        pair = _canonical_expression(tuple(self.pair_terms))
        previous = _canonical_expression(tuple(self.previous_shifted_terms))
        residual = _canonical_expression(linear + pair + _negated(previous))
        if residual:
            rendered = ", ".join(f"{term.atom}:{term.coefficient}" for term in residual)
            raise ValueError(f"retained recurrence is not an exact zero identity: {rendered}")

        object.__setattr__(self, "order", order)
        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "coefficient_state_id", coefficient_state_id)
        object.__setattr__(self, "identity_source", identity_source)
        object.__setattr__(self, "theorem_name", theorem_name)
        object.__setattr__(self, "dependencies", dependencies)
        object.__setattr__(self, "linear_terms", linear)
        object.__setattr__(self, "pair_terms", pair)
        object.__setattr__(self, "previous_shifted_terms", previous)

    @property
    def exact_zero(self) -> bool:
        return True

    @property
    def residual_terms(self) -> tuple[FormalAtomTerm, ...]:
        return ()


@dataclass(frozen=True)
class FiniteRetainedCancellationCertificate:
    """Contiguous exact retained-cancellation prefix ``n=0,...,N``."""

    hierarchy_id: str
    source_revision: str
    coefficient_state_id: str
    max_order: int
    identities: tuple[ExactRetainedRecurrenceIdentity, ...]

    def __post_init__(self) -> None:
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        coefficient_state_id = _nonempty_text(self.coefficient_state_id, "coefficient_state_id")
        max_order = _nonnegative_int(self.max_order, "max_order")
        identities = tuple(self.identities)
        if not all(isinstance(identity, ExactRetainedRecurrenceIdentity) for identity in identities):
            raise TypeError("identities must contain ExactRetainedRecurrenceIdentity values")
        orders = [identity.order for identity in identities]
        if len(orders) != len(set(orders)):
            raise ValueError("retained cancellation certificate contains duplicate orders")
        expected = list(range(max_order + 1))
        if sorted(orders) != expected:
            raise ValueError(
                f"retained cancellation certificate must cover the contiguous prefix {expected}"
            )
        for identity in identities:
            if identity.hierarchy_id != hierarchy_id:
                raise ValueError("retained identity hierarchy_id does not match its certificate")
            if identity.source_revision != source_revision:
                raise ValueError("retained identity source_revision does not match its certificate")
            if identity.coefficient_state_id != coefficient_state_id:
                raise ValueError("retained identity coefficient_state_id does not match its certificate")
            if not identity.exact_zero:
                raise ValueError("every retained identity must be exact zero")
        ordered = tuple(sorted(identities, key=lambda identity: identity.order))
        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "coefficient_state_id", coefficient_state_id)
        object.__setattr__(self, "max_order", max_order)
        object.__setattr__(self, "identities", ordered)

    @property
    def paper_exact(self) -> bool:
        return False

    def identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
        order = _nonnegative_int(order, "order")
        if order > self.max_order:
            raise ValueError("requested retained identity lies outside the certified prefix")
        identity = self.identities[order]
        if identity.order != order:
            raise RuntimeError("internal retained-cancellation ordering invariant was violated")
        return identity


@dataclass(frozen=True)
class CertifiedFullResidualTailMajorant:
    """First-omitted-order bound admitted only after exact retained cancellation."""

    cancellations: FiniteRetainedCancellationCertificate
    tail: SlowTailMajorant

    def __post_init__(self) -> None:
        if not isinstance(self.cancellations, FiniteRetainedCancellationCertificate):
            raise TypeError("cancellations must be a FiniteRetainedCancellationCertificate")
        if not isinstance(self.tail, SlowTailMajorant):
            raise TypeError("tail must be a SlowTailMajorant")
        if self.tail.order != self.cancellations.max_order:
            raise ValueError("tail order does not match retained cancellation prefix")

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def first_omitted_order(self) -> int:
        return self.tail.order + 1


def certify_full_residual_tail_majorant(
    cancellations: FiniteRetainedCancellationCertificate,
    q: float,
    h: float,
    base_power: float,
    pair: Sequence[Sequence[float]],
    shifted: Sequence[float],
) -> CertifiedFullResidualTailMajorant:
    """Expose the omitted-tail majorant only behind an exact cancellation prefix.

    The pair/shifted coefficients remain numerical inputs to the finite tail
    majorant, so this is not yet an all-order or paper-exact residual proof.
    Its purpose is narrower: callers can no longer promote the omitted-tail
    majorant to a full-residual statement merely because floating recurrence
    values happen to be small.
    """

    if not isinstance(cancellations, FiniteRetainedCancellationCertificate):
        raise TypeError("cancellations must be a FiniteRetainedCancellationCertificate")
    tail = omitted_slow_tail_majorant(
        q,
        h,
        base_power,
        cancellations.max_order,
        pair,
        shifted,
    )
    return CertifiedFullResidualTailMajorant(cancellations=cancellations, tail=tail)
