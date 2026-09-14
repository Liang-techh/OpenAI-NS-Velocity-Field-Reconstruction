"""Exact arithmetic certificate for the Section 8 five-row mean repair algebra.

This module replays only the finite-dimensional identities pinned in
``NavierStokes/FiveRowRank.lean`` and ``NavierStokes/MeanRankUpdate.lean``.
It deliberately does not manufacture the noncomputable localized smooth bumps,
materialize ``CorrectionState.debt``, or identify a caller-supplied defect with
the actual oscillatory-wave defect.

The certificate uses :class:`fractions.Fraction` throughout.  In particular,
zero auxiliary moments and the three debt-cancelling rows are checked as exact
equalities, never by a numerical tolerance.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Sequence


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FIVE_ROW_FILE = "NavierStokes/FiveRowRank.lean"
PINNED_SCALE_FILE = "NavierStokes/MeanRankUpdate.lean"
PINNED_FIVE_ROWS_SYMBOL = "NavierStokes.FiveRowRank.FiveRows"
PINNED_NORMALIZE_DEBT_SYMBOL = "NavierStokes.MeanRankUpdate.normalizeDebt"
PINNED_PHYSICAL_FIVE_ROWS_SYMBOL = "NavierStokes.MeanRankUpdate.physical_five_rows"


def _exact(value: object, name: str) -> Fraction:
    """Embed integers/Fractions exactly and reject approximate scalar inputs."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an int or Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(
        f"{name} must be an int or Fraction; float/Decimal approximations are rejected"
    )


def _vector(values: Sequence[object], length: int, name: str) -> tuple[Fraction, ...]:
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} entries")
    return tuple(_exact(value, f"{name}[{index}]") for index, value in enumerate(values))


@dataclass(frozen=True)
class ExactFiveRowAlgebraCertificate:
    """Exact normalized/physical algebra for one Section 8 debt vector.

    ``debt`` is the physical ``(P, J_theta, J_z)`` debt.  The formal localized
    repair first normalizes this debt by the physical length/velocity scales,
    imposes one zero angular auxiliary moment and one zero axial auxiliary
    moment, and uses the remaining three moments to cancel the three debt rows.

    This class certifies exactly that algebraic implication.  It does *not*
    certify existence/evaluation of the localized bump inverse itself.
    """

    lam: Fraction | int
    C: Fraction | int
    ell: Fraction | int
    U: Fraction | int
    debt: tuple[Fraction | int, Fraction | int, Fraction | int]

    def __post_init__(self) -> None:
        lam = _exact(self.lam, "lam")
        C = _exact(self.C, "C")
        ell = _exact(self.ell, "ell")
        U = _exact(self.U, "U")
        debt = _vector(self.debt, 3, "debt")
        if lam <= 0:
            raise ValueError("lam must be positive")
        if C == 0:
            raise ValueError("C must be nonzero")
        if ell <= 0:
            raise ValueError("ell must be positive")
        if U == 0:
            raise ValueError("U must be nonzero")
        object.__setattr__(self, "lam", lam)
        object.__setattr__(self, "C", C)
        object.__setattr__(self, "ell", ell)
        object.__setattr__(self, "U", U)
        object.__setattr__(self, "debt", debt)
        if len(set(self.angular_powers)) != 3:
            raise ValueError("positive lam must give three distinct angular powers")
        if len(set(self.axial_powers)) != 2:
            raise ValueError("positive lam must give two distinct axial powers")
        if not self.five_row_identity_holds:
            raise ValueError("pinned five-row algebra identity failed")

    @property
    def angular_powers(self) -> tuple[Fraction, Fraction, Fraction]:
        lam = self.lam
        return Fraction(2), -Fraction(2) - 2 * lam, -2 * lam

    @property
    def axial_powers(self) -> tuple[Fraction, Fraction]:
        lam = self.lam
        return Fraction(1), Fraction(1) - 2 * lam

    @property
    def normalized_debt(self) -> tuple[Fraction, Fraction, Fraction]:
        d0, d1, d2 = self.debt
        return (
            d0 / self.U**2,
            d1 / (self.ell**3 * self.U**2),
            d2 / (self.ell**2 * self.U**2),
        )

    @property
    def angular_moment_targets(self) -> tuple[Fraction, Fraction, Fraction]:
        d0, _, d2 = self.normalized_debt
        return Fraction(0), -d0 / (2 * self.C), d2 / self.C

    @property
    def axial_moment_targets(self) -> tuple[Fraction, Fraction]:
        _, d1, _ = self.normalized_debt
        return Fraction(0), -d1 / self.C

    @property
    def zero_auxiliary_average_identity_holds(self) -> bool:
        return self.angular_moment_targets[0] == 0 and self.axial_moment_targets[0] == 0

    @property
    def physical_row_targets(self) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
        d0, d1, d2 = self.debt
        return Fraction(0), Fraction(0), -d0, -d1, -d2

    def physical_rows_from_normalized_moments(
        self,
        angular: Sequence[object],
        axial: Sequence[object],
    ) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
        """Replay the five exact physical scaling rows from normalized moments."""
        am = _vector(angular, 3, "angular")
        ax = _vector(axial, 2, "axial")
        return (
            self.ell**3 * self.U * am[0],
            self.ell**2 * self.U * ax[0],
            self.U**2 * (2 * self.C * am[1]),
            self.ell**3 * self.U**2 * (self.C * ax[1]),
            self.ell**2 * self.U**2 * (-self.C * am[2]),
        )

    @property
    def certified_physical_rows(self) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
        return self.physical_rows_from_normalized_moments(
            self.angular_moment_targets,
            self.axial_moment_targets,
        )

    @property
    def five_row_identity_holds(self) -> bool:
        return (
            self.zero_auxiliary_average_identity_holds
            and self.certified_physical_rows == self.physical_row_targets
        )

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def exact_five_row_algebra_machine_checked(self) -> bool:
        return True

    @property
    def localized_moment_inverse_materialized(self) -> bool:
        return False

    @property
    def actual_wave_defect_consumed(self) -> bool:
        return False

    @property
    def actual_compact_correction_field_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
