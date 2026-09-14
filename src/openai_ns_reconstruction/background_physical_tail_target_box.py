"""Finite all-jet / integer-power target boxes for Section 5 physical tails.

The pinned all-order SlowBorel argument has the quantifier shape that every
finite physical derivative order and every requested algebraic decay power can
be reached by a sufficiently late prefix of one recursively chosen diagonal
schedule. The preceding cofinal-ladder layer materializes finitely many
canonical diagonal requests ``level k := (M=k, P=k)`` and proves those levels
are literal extensions of one provider-owned coefficient/schedule family.

This module packages the next finite quantifier shape without promoting it to
an infinite theorem: for one rectangle

    0 <= m <= M,   1 <= N <= P,

one single canonical frontier level ``K=max(M,P,1)`` must dominate every cell.
Each cell is tied to the exact derivative row of that same frontier certificate
and records its least canonical level ``max(m,N,1)``. No new C[j,m] rows,
support witnesses, recurrence identities, cutoff scales, physical-chart
constants, or residual values are caller supplied here.

A finite target box is still finite. It does not prove provider totality, an
infinite diagonal schedule, an actual PDE residual estimate, all-jets flatness,
or super-algebraic convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .background_physical_tail_cofinal_ladder import (
    FiniteCofinalPhysicalTailLadderCertificate,
    canonical_level_for_physical_target,
)
from .background_physical_tail_prefactor import (
    FinitePhysicalTailPrefactorCertificate,
    PhysicalTailPrefactorRow,
)


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _positive_int(value: int, name: str) -> int:
    result = _nonnegative_int(value, name)
    if result == 0:
        raise ValueError(f"{name} must be positive")
    return result


@dataclass(frozen=True)
class FinitePhysicalTailTargetCell:
    """One exact ``(m,N)`` request witnessed by the common frontier prefix."""

    derivative_order: int
    requested_integer_power: int
    minimal_canonical_level: int
    frontier_level: int
    frontier_certificate: FinitePhysicalTailPrefactorCertificate
    prefactor_row: PhysicalTailPrefactorRow

    def __post_init__(self) -> None:
        m = _nonnegative_int(self.derivative_order, "derivative_order")
        N = _positive_int(self.requested_integer_power, "requested_integer_power")
        minimal = _positive_int(self.minimal_canonical_level, "minimal_canonical_level")
        frontier = _positive_int(self.frontier_level, "frontier_level")
        if not isinstance(self.frontier_certificate, FinitePhysicalTailPrefactorCertificate):
            raise TypeError(
                "frontier_certificate must be a FinitePhysicalTailPrefactorCertificate"
            )
        if not isinstance(self.prefactor_row, PhysicalTailPrefactorRow):
            raise TypeError("prefactor_row must be a PhysicalTailPrefactorRow")

        expected_minimal = canonical_level_for_physical_target(m, Fraction(N, 1))
        if minimal != expected_minimal:
            raise ValueError("minimal_canonical_level is not max(m,N,1)")
        if frontier < minimal:
            raise ValueError("common frontier level does not dominate this target cell")

        physical = self.frontier_certificate.first_omitted.uncut.physical
        if physical.max_derivative_order != frontier:
            raise ValueError("frontier certificate derivative budget is not the frontier level")
        if physical.target_physical_power != Fraction(frontier, 1):
            raise ValueError("frontier certificate target power is not the frontier level")

        rows = self.frontier_certificate.rows
        if m >= len(rows):
            raise ValueError("frontier certificate does not contain the requested derivative row")
        if self.prefactor_row != rows[m]:
            raise ValueError("target cell is not bound to the frontier certificate derivative row")
        if self.prefactor_row.derivative_order != m:
            raise RuntimeError("frontier prefactor row derivative order is misaligned")
        if self.prefactor_row.physical_power < N:
            raise ValueError("frontier physical tail power does not dominate the requested power")
        if self.prefactor_row.dyadic_prefactor_exact != (
            self.frontier_certificate.first_omitted.dyadic_prefactor_exact
        ):
            raise ValueError("target cell lost the frontier dyadic prefactor")

    @property
    def requested_power_exact(self) -> Fraction:
        return Fraction(self.requested_integer_power, 1)

    @property
    def physical_power_slack(self) -> Fraction:
        return self.prefactor_row.physical_power - self.requested_power_exact


@dataclass(frozen=True)
class FinitePhysicalTailTargetBoxCertificate:
    """One common finite prefix for all ``m<=M`` and integer powers ``N<=P``."""

    ladder: FiniteCofinalPhysicalTailLadderCertificate
    max_derivative_order: int
    max_integer_power: int
    frontier_level: int
    frontier_certificate: FinitePhysicalTailPrefactorCertificate
    cells: tuple[FinitePhysicalTailTargetCell, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.ladder, FiniteCofinalPhysicalTailLadderCertificate):
            raise TypeError("ladder must be a FiniteCofinalPhysicalTailLadderCertificate")
        M = _nonnegative_int(self.max_derivative_order, "max_derivative_order")
        P = _positive_int(self.max_integer_power, "max_integer_power")
        frontier = _positive_int(self.frontier_level, "frontier_level")
        expected_frontier = max(M, P, 1)
        if frontier != expected_frontier:
            raise ValueError("frontier_level must equal max(M,P,1)")
        if frontier > self.ladder.max_level:
            raise ValueError("finite target box exceeds the materialized cofinal ladder")
        if not isinstance(self.frontier_certificate, FinitePhysicalTailPrefactorCertificate):
            raise TypeError(
                "frontier_certificate must be a FinitePhysicalTailPrefactorCertificate"
            )

        # Object identity is intentional: the box must consume the already
        # materialized canonical ladder level, not an equal-looking side input.
        if self.frontier_certificate is not self.ladder.levels[frontier - 1]:
            raise ValueError("frontier certificate is not the ladder's materialized frontier level")

        expected_pairs = tuple((m, N) for m in range(M + 1) for N in range(1, P + 1))
        if len(self.cells) != len(expected_pairs):
            raise ValueError("finite target box has the wrong number of target cells")
        actual_pairs = tuple(
            (cell.derivative_order, cell.requested_integer_power) for cell in self.cells
        )
        if actual_pairs != expected_pairs:
            raise ValueError("finite target box must cover every (m,N) cell exactly once")

        for cell in self.cells:
            if not isinstance(cell, FinitePhysicalTailTargetCell):
                raise TypeError("cells must contain FinitePhysicalTailTargetCell values")
            if cell.frontier_level != frontier:
                raise ValueError("target cell changed the common frontier level")
            if cell.frontier_certificate is not self.frontier_certificate:
                raise ValueError("target cell is cross-wired to a different frontier certificate")
            expected_minimal = canonical_level_for_physical_target(
                cell.derivative_order, Fraction(cell.requested_integer_power, 1)
            )
            if cell.minimal_canonical_level != expected_minimal:
                raise ValueError("target cell minimal canonical level drifted")

    @property
    def target_count(self) -> int:
        return len(self.cells)

    def cell(
        self, derivative_order: int, requested_integer_power: int
    ) -> FinitePhysicalTailTargetCell:
        m = _nonnegative_int(derivative_order, "derivative_order")
        N = _positive_int(requested_integer_power, "requested_integer_power")
        if m > self.max_derivative_order or N > self.max_integer_power:
            raise ValueError("requested cell lies outside the certified finite target box")
        index = m * self.max_integer_power + (N - 1)
        result = self.cells[index]
        if (result.derivative_order, result.requested_integer_power) != (m, N):
            raise RuntimeError("finite target box cell indexing drifted")
        return result

    @property
    def one_common_frontier_prefix_verified(self) -> bool:
        return True

    @property
    def finite_all_jet_integer_power_box_verified(self) -> bool:
        return True

    @property
    def finite_target_box_only(self) -> bool:
        return True

    @property
    def provider_totality_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def actual_pde_residual_verified(self) -> bool:
        return False

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def all_jets_flat(self) -> bool:
        return False

    @property
    def super_algebraic(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False


def certify_finite_physical_tail_target_box(
    ladder: FiniteCofinalPhysicalTailLadderCertificate,
    *,
    max_derivative_order: int,
    max_integer_power: int,
) -> FinitePhysicalTailTargetBoxCertificate:
    """Certify one common finite prefix for a rectangular jet/power target set."""

    if not isinstance(ladder, FiniteCofinalPhysicalTailLadderCertificate):
        raise TypeError("ladder must be a FiniteCofinalPhysicalTailLadderCertificate")
    M = _nonnegative_int(max_derivative_order, "max_derivative_order")
    P = _positive_int(max_integer_power, "max_integer_power")
    frontier = max(M, P, 1)
    if frontier > ladder.max_level:
        raise ValueError("finite target box exceeds the materialized cofinal ladder")

    frontier_certificate = ladder.levels[frontier - 1]
    rows = frontier_certificate.rows
    cells = tuple(
        FinitePhysicalTailTargetCell(
            derivative_order=m,
            requested_integer_power=N,
            minimal_canonical_level=canonical_level_for_physical_target(m, Fraction(N, 1)),
            frontier_level=frontier,
            frontier_certificate=frontier_certificate,
            prefactor_row=rows[m],
        )
        for m in range(M + 1)
        for N in range(1, P + 1)
    )
    return FinitePhysicalTailTargetBoxCertificate(
        ladder=ladder,
        max_derivative_order=M,
        max_integer_power=P,
        frontier_level=frontier,
        frontier_certificate=frontier_certificate,
        cells=cells,
    )
