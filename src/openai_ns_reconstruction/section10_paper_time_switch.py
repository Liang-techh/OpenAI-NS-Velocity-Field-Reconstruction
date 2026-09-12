"""Pinned Section 10 time-switch provenance and exact theorem geometry.

This module intentionally does *not* implement a numerical surrogate for the
paper's time cutoff.  The pinned formal source defines ``timeSwitch`` using a
Mathlib ``ContDiffBump`` and proves exact plateau/support facts.  We record only
those theorem-facing facts so downstream code cannot silently substitute an
arbitrary ``TimeWindowCutoff`` and call it the paper localization.

The source also proves that ``timeSwitch`` is identically one after ``3/4``
and that every positive-order derivative vanishes for ``t > 3/4``.  In
particular the switch is inert in a neighborhood of ``t = 1``.  This is *not*
the missing smooth extension of the Section 9 field through ``t = 1``; that
still requires the actual one-sided field/residual construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Optional


PINNED_LEAN_REPOSITORY = "openai/NavierStokesAndEuler"
PINNED_LEAN_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_SOURCE_FILE = "NavierStokes/SmoothCutoffs.lean"
PINNED_TIME_SWITCH_DEFINITION = "NavierStokes.SmoothCutoffs.timeSwitch"
PINNED_TIME_SWITCH_CONTDIFF = "NavierStokes.SmoothCutoffs.timeSwitch_contDiff"
PINNED_ZERO_THEOREM = "NavierStokes.SmoothCutoffs.timeSwitch_zero_of_abs_le"
PINNED_ONE_THEOREM = (
    "NavierStokes.SmoothCutoffs.timeSwitch_one_of_three_quarters_le"
)
PINNED_EVENTUALLY_ONE_THEOREM = (
    "NavierStokes.SmoothCutoffs.timeSwitch_eventually_one"
)
PINNED_LATE_DERIVATIVE_THEOREM = (
    "NavierStokes.SmoothCutoffs.timeSwitch_iteratedDeriv_late"
)
PINNED_DERIVATIVE_SUPPORT_THEOREM = (
    "NavierStokes.SmoothCutoffs.timeSwitch_iteratedDeriv_support_nonneg"
)

ZERO_PLATEAU_ABS_RADIUS = Fraction(3, 8)
ONE_PLATEAU_START = Fraction(3, 4)


def _as_exact_rational(value: object, name: str) -> Fraction:
    """Accept only exact integer/Fraction input for theorem-boundary queries."""

    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact rational, not bool")
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, Fraction):
        return value
    raise TypeError(
        f"{name} must be int or fractions.Fraction; float input is rejected "
        "at exact theorem boundaries"
    )


@dataclass(frozen=True)
class Section10PaperTimeSwitchSource:
    """Fail-closed binding to the paper's pinned formal ``timeSwitch`` source.

    The caller supplies provenance metadata, but every field is checked against
    the pinned revision and exact theorem geometry.  This prevents a generic
    user-selected time window from being admitted as the Section 10 switch.
    """

    lean_repository: str
    lean_commit: str
    source_file: str
    definition_name: str
    contdiff_theorem: str
    zero_theorem: str
    one_theorem: str
    eventually_one_theorem: str
    late_derivative_theorem: str
    derivative_support_theorem: str
    zero_plateau_abs_radius: Fraction
    one_plateau_start: Fraction

    def __post_init__(self) -> None:
        expected = {
            "lean_repository": PINNED_LEAN_REPOSITORY,
            "lean_commit": PINNED_LEAN_COMMIT,
            "source_file": PINNED_SOURCE_FILE,
            "definition_name": PINNED_TIME_SWITCH_DEFINITION,
            "contdiff_theorem": PINNED_TIME_SWITCH_CONTDIFF,
            "zero_theorem": PINNED_ZERO_THEOREM,
            "one_theorem": PINNED_ONE_THEOREM,
            "eventually_one_theorem": PINNED_EVENTUALLY_ONE_THEOREM,
            "late_derivative_theorem": PINNED_LATE_DERIVATIVE_THEOREM,
            "derivative_support_theorem": PINNED_DERIVATIVE_SUPPORT_THEOREM,
        }
        for name, wanted in expected.items():
            if getattr(self, name) != wanted:
                raise ValueError(f"{name} must match the pinned SmoothCutoffs source")

        if not isinstance(self.zero_plateau_abs_radius, Fraction):
            raise TypeError("zero_plateau_abs_radius must be fractions.Fraction")
        if not isinstance(self.one_plateau_start, Fraction):
            raise TypeError("one_plateau_start must be fractions.Fraction")
        if self.zero_plateau_abs_radius != ZERO_PLATEAU_ABS_RADIUS:
            raise ValueError("zero plateau must be the theorem value 3/8")
        if self.one_plateau_start != ONE_PLATEAU_START:
            raise ValueError("one plateau must start at the theorem value 3/4")

    @classmethod
    def pinned(cls) -> "Section10PaperTimeSwitchSource":
        return cls(
            lean_repository=PINNED_LEAN_REPOSITORY,
            lean_commit=PINNED_LEAN_COMMIT,
            source_file=PINNED_SOURCE_FILE,
            definition_name=PINNED_TIME_SWITCH_DEFINITION,
            contdiff_theorem=PINNED_TIME_SWITCH_CONTDIFF,
            zero_theorem=PINNED_ZERO_THEOREM,
            one_theorem=PINNED_ONE_THEOREM,
            eventually_one_theorem=PINNED_EVENTUALLY_ONE_THEOREM,
            late_derivative_theorem=PINNED_LATE_DERIVATIVE_THEOREM,
            derivative_support_theorem=PINNED_DERIVATIVE_SUPPORT_THEOREM,
            zero_plateau_abs_radius=ZERO_PLATEAU_ABS_RADIUS,
            one_plateau_start=ONE_PLATEAU_START,
        )

    @property
    def source_key(self) -> tuple[str, str, str]:
        return self.lean_repository, self.lean_commit, self.source_file

    def certified_value(self, t: int | Fraction) -> Optional[Fraction]:
        """Return only values proved by the pinned plateau theorems.

        ``None`` means that the formal source does not force a plateau value at
        that time.  The transition is deliberately not numerically evaluated.
        """

        q = _as_exact_rational(t, "t")
        if abs(q) <= self.zero_plateau_abs_radius:
            return Fraction(0, 1)
        if q >= self.one_plateau_start:
            return Fraction(1, 1)
        return None

    def positive_derivatives_certified_zero_late(self, t: int | Fraction) -> bool:
        """Whether the pinned all-order late-derivative theorem applies.

        The Lean theorem uses the strict hypothesis ``3/4 < t``.  We preserve
        that strict boundary instead of inferring an endpoint derivative fact.
        """

        q = _as_exact_rational(t, "t")
        return q > self.one_plateau_start

    def positive_derivative_support_may_be_nonzero_on_nonnegative_axis(
        self, t: int | Fraction
    ) -> bool:
        """Necessary support collar from the pinned theorem, for ``t >= 0``.

        This answers only whether a positive derivative *may* be nonzero from
        the theorem support enclosure; it does not claim non-vanishing.
        """

        q = _as_exact_rational(t, "t")
        if q < 0:
            raise ValueError("the derivative-support theorem is restricted to t >= 0")
        return self.zero_plateau_abs_radius <= q <= self.one_plateau_start

    @property
    def time_switch_contdiff_theorem_bound(self) -> bool:
        return self.contdiff_theorem == PINNED_TIME_SWITCH_CONTDIFF

    @property
    def endpoint_t1_switch_value_certified_one(self) -> bool:
        return self.certified_value(Fraction(1, 1)) == 1

    @property
    def endpoint_t1_switch_positive_derivatives_certified_zero(self) -> bool:
        return self.positive_derivatives_certified_zero_late(Fraction(1, 1))

    @property
    def actual_mathlib_bump_numerically_evaluated(self) -> bool:
        return False

    @property
    def lean_theorems_machine_replayed_in_python(self) -> bool:
        return False

    @property
    def section9_field_smooth_extension_through_t1_constructed(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
