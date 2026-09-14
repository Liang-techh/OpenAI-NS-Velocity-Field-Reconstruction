"""Definition-exact support geometry for the Section 8 five-row repair.

The paper/formal construction does not choose arbitrary compact cutoffs.  In
``NavierStokes/FiveRowRank.lean`` the three angular and two axial correction
cells are fixed affine subintervals of a positive patch ``(a,b)``.  The actual
smooth repairs from ``LocalizedMomentRepair.lean`` then live in the middle half
of each cell.

The physical endpoints ``a`` and ``b`` can be noncomputable real expressions,
so this module never approximates them.  Instead it works in the exact affine
coordinate ``R = (1-alpha) * a + alpha * b``.  All support/order statements
therefore reduce to rational inequalities in ``alpha`` and can be replayed
without floating point, sampling, quadrature, or a caller-supplied cutoff.

This certificate closes only the support-geometry part of Section 8.  It does
not materialize the Lean bump, its moment matrix, the five-row correction, or
an oscillatory-wave defect, and it is not a paper-exact velocity certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Iterable


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FIVE_ROW_FILE = "NavierStokes/FiveRowRank.lean"
PINNED_REPAIR_FILE = "NavierStokes/LocalizedMomentRepair.lean"
PINNED_SYMBOLS = (
    "NavierStokes.FiveRowRank.cellStep",
    "NavierStokes.FiveRowRank.cellLower",
    "NavierStokes.FiveRowRank.cellUpper",
    "NavierStokes.FiveRowRank.cell_separated",
    "NavierStokes.FiveRowRank.cell_union_subset",
    "NavierStokes.LocalizedMomentRepair.innerLower",
    "NavierStokes.LocalizedMomentRepair.innerUpper",
    "NavierStokes.LocalizedMomentRepair.repair_tsupport_subset_open",
)
ANGULAR_CELL_COUNT = 3
AXIAL_CELL_COUNT = 2


def _count(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("count must be a positive integer")
    return int(value)


def _exact_fraction(value: object, name: str) -> Fraction:
    """Admit only exact integer/rational affine coordinates.

    In particular, binary floats are rejected even when they happen to display
    as a simple rational.  A support certificate must not depend on a rounded
    representation of the pinned affine formulas.
    """
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(f"{name} must be an exact integer/Fraction; floats are rejected")


def pinned_cell_coordinates(count: int) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return ``cellLower/cellUpper`` in exact affine patch coordinates."""
    n = _count(count)
    denominator = 2 * n + 1
    return tuple(
        (
            Fraction(2 * j + 1, denominator),
            Fraction(2 * j + 2, denominator),
        )
        for j in range(n)
    )


def pinned_inner_support(
    cell: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    """Return ``innerLower/innerUpper`` for one cell in affine coordinates."""
    if len(cell) != 2:
        raise ValueError("cell must contain exactly two affine coordinates")
    left = _exact_fraction(cell[0], "cell lower")
    right = _exact_fraction(cell[1], "cell upper")
    if not left < right:
        raise ValueError("cell lower must be strictly below cell upper")
    return (3 * left + right) / 4, (left + 3 * right) / 4


def _normalize_pairs(
    pairs: Iterable[tuple[object, object]], name: str
) -> tuple[tuple[Fraction, Fraction], ...]:
    normalized = []
    for i, pair in enumerate(pairs):
        if len(pair) != 2:
            raise ValueError(f"{name}[{i}] must contain exactly two coordinates")
        normalized.append(
            (
                _exact_fraction(pair[0], f"{name}[{i}].lower"),
                _exact_fraction(pair[1], f"{name}[{i}].upper"),
            )
        )
    return tuple(normalized)


@dataclass(frozen=True)
class AffineSupportFamily:
    """Exact replay of one pinned ``FiveRowRank`` cell/support family."""

    count: int
    cells: tuple[tuple[Fraction, Fraction], ...]
    supports: tuple[tuple[Fraction, Fraction], ...]

    def __post_init__(self) -> None:
        n = _count(self.count)
        object.__setattr__(self, "count", n)
        cells = _normalize_pairs(self.cells, "cells")
        supports = _normalize_pairs(self.supports, "supports")
        object.__setattr__(self, "cells", cells)
        object.__setattr__(self, "supports", supports)

        expected_cells = pinned_cell_coordinates(n)
        if cells != expected_cells:
            raise ValueError("cells do not equal the pinned FiveRowRank cellLower/cellUpper geometry")
        expected_supports = tuple(pinned_inner_support(cell) for cell in expected_cells)
        if supports != expected_supports:
            raise ValueError(
                "supports do not equal the pinned LocalizedMomentRepair middle-half geometry"
            )

        for i, ((left, right), (support_left, support_right)) in enumerate(
            zip(cells, supports)
        ):
            if not Fraction(0) < left < support_left < support_right < right < Fraction(1):
                raise ValueError(f"cell/support {i} is not strictly inside the normalized patch")
        for i in range(n - 1):
            if not cells[i][1] < cells[i + 1][0]:
                raise ValueError(f"cells {i} and {i + 1} are not separated by a positive gap")

    @classmethod
    def pinned(cls, count: int) -> "AffineSupportFamily":
        cells = pinned_cell_coordinates(count)
        return cls(
            count=count,
            cells=cells,
            supports=tuple(pinned_inner_support(cell) for cell in cells),
        )

    @property
    def exact_affine_geometry_verified(self) -> bool:
        return True

    @property
    def pairwise_separated_verified(self) -> bool:
        return True

    @property
    def support_strictly_inside_cells_verified(self) -> bool:
        return True

    @property
    def cells_strictly_inside_patch_verified(self) -> bool:
        return True


@dataclass(frozen=True)
class Section8CompactSupportGeometry:
    """Pinned three-angular/two-axial Section 8 support certificate."""

    angular: AffineSupportFamily
    axial: AffineSupportFamily
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    dependency_files: tuple[str, str] = (PINNED_FIVE_ROW_FILE, PINNED_REPAIR_FILE)
    dependency_symbols: tuple[str, ...] = PINNED_SYMBOLS

    def __post_init__(self) -> None:
        if not isinstance(self.angular, AffineSupportFamily):
            raise TypeError("angular must be an AffineSupportFamily")
        if not isinstance(self.axial, AffineSupportFamily):
            raise TypeError("axial must be an AffineSupportFamily")
        if self.angular.count != ANGULAR_CELL_COUNT:
            raise ValueError("Section 8 angular family must use exactly three cells")
        if self.axial.count != AXIAL_CELL_COUNT:
            raise ValueError("Section 8 axial family must use exactly two cells")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository does not match the pinned source")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit does not match the pinned source")
        if self.dependency_files != (PINNED_FIVE_ROW_FILE, PINNED_REPAIR_FILE):
            raise ValueError("dependency_files do not match the pinned support chain")
        if self.dependency_symbols != PINNED_SYMBOLS:
            raise ValueError("dependency_symbols do not match the pinned support chain")

    @classmethod
    def pinned(cls) -> "Section8CompactSupportGeometry":
        return cls(
            angular=AffineSupportFamily.pinned(ANGULAR_CELL_COUNT),
            axial=AffineSupportFamily.pinned(AXIAL_CELL_COUNT),
        )

    @property
    def status(self) -> str:
        return "definition-exact-support-geometry"

    @property
    def valid_for_any_strict_positive_patch(self) -> bool:
        """The affine proof applies after supplying the paper hypotheses ``0<a<b``."""
        return True

    @property
    def actual_patch_endpoints_materialized(self) -> bool:
        return False

    @property
    def actual_bump_values_materialized(self) -> bool:
        return False

    @property
    def exact_moment_identities_machine_replayed(self) -> bool:
        return False

    @property
    def five_row_solve_materialized(self) -> bool:
        return False

    @property
    def actual_wave_defect_consumed(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
