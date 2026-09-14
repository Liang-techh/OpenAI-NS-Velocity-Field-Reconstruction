"""Symbolic physical-chart prefactor ledger for finite Section 5 tail requests.

The pinned proof of ``SlowBorelBase.exists_physical_uncut_tail`` does not need a
numerically evaluated physical-coordinate constant.  For a fixed derivative
budget ``m`` it first obtains positive constants ``C_i`` from
``physicalChart_jet_bound``, defines

    D = 1 + sum_{i=0}^m C_i,

and then the chain rule produces the finite prefactor

    m! * 2^-J * D^m

in front of the requested physical q-power.

This module records that exact dependency graph without inventing numerical
values for the existential ``C_i``.  The support interval, truncation order,
dyadic factor, target power, and first-omitted exponent all come from the
provider-owned certificate built upstream.  Callers cannot supply ``D``,
``C_i``, another support interval, or another cutoff schedule.

The result is theorem-shape bookkeeping only: Python does not replay Lean,
materialize the existential real constants, evaluate the PDE residual, prove a
total all-order coefficient hierarchy, or establish all-jets flatness.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial

from .background_first_omitted_physical_tail import (
    FiniteFirstOmittedPhysicalTailCertificate,
    certify_finite_first_omitted_physical_tail,
)
from .background_target_driven_exact_majorants import TargetDrivenExactMajorantProvider


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_CHART_JET_THEOREM = "NavierStokes.SlowBorelBase.physicalChart_jet_bound"
PINNED_FINITE_BOUND_THEOREM = "NavierStokes.SlowBorelBase.physicalChart_finite_bound"
PINNED_COMPOSITION_THEOREM = "NavierStokes.SlowBorelBase.physical_composition_bound"
PINNED_PHYSICAL_TAIL_THEOREM = "NavierStokes.SlowBorelBase.exists_physical_uncut_tail"


@dataclass(frozen=True)
class PhysicalChartExistentialJetConstant:
    """One pinned existential ``C_i>0`` used to assemble the finite ``D``."""

    derivative_order: int
    symbol: str
    formal_revision: str = PINNED_FORMAL_REVISION
    theorem_name: str = PINNED_CHART_JET_THEOREM

    def __post_init__(self) -> None:
        if isinstance(self.derivative_order, bool) or not isinstance(self.derivative_order, int):
            raise TypeError("derivative_order must be an integer")
        if self.derivative_order < 0:
            raise ValueError("derivative_order must be nonnegative")
        expected = f"C_phys_{self.derivative_order}"
        if self.symbol != expected:
            raise ValueError(f"physical-chart existential symbol must be {expected!r}")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.theorem_name != PINNED_CHART_JET_THEOREM:
            raise ValueError("theorem_name does not match physicalChart_jet_bound")


@dataclass(frozen=True)
class PhysicalChartFiniteBoundShape:
    """Pinned finite ``D=1+sum C_i`` shape on the hierarchy-owned X interval."""

    max_derivative_order: int
    h_exact: Fraction
    x_lower: Fraction
    x_upper: Fraction
    jet_constants: tuple[PhysicalChartExistentialJetConstant, ...]
    formal_revision: str = PINNED_FORMAL_REVISION
    theorem_name: str = PINNED_FINITE_BOUND_THEOREM

    def __post_init__(self) -> None:
        if isinstance(self.max_derivative_order, bool) or not isinstance(
            self.max_derivative_order, int
        ):
            raise TypeError("max_derivative_order must be an integer")
        if self.max_derivative_order < 0:
            raise ValueError("max_derivative_order must be nonnegative")
        if not isinstance(self.h_exact, Fraction):
            raise TypeError("h_exact must be an exact Fraction")
        if not Fraction(0, 1) < self.h_exact < Fraction(1, 2):
            raise ValueError("pinned physical-chart theorem requires 0 < h < 1/2")
        if not isinstance(self.x_lower, Fraction) or not isinstance(self.x_upper, Fraction):
            raise TypeError("physical-chart X interval must use exact Fraction endpoints")
        if self.x_lower > self.x_upper:
            raise ValueError("physical-chart X interval is reversed")
        constants = tuple(self.jet_constants)
        if not all(isinstance(c, PhysicalChartExistentialJetConstant) for c in constants):
            raise TypeError("jet_constants must contain PhysicalChartExistentialJetConstant values")
        expected_orders = tuple(range(self.max_derivative_order + 1))
        if tuple(c.derivative_order for c in constants) != expected_orders:
            raise ValueError("physical-chart constants must cover exactly 0..M")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.theorem_name != PINNED_FINITE_BOUND_THEOREM:
            raise ValueError("theorem_name does not match physicalChart_finite_bound")

    @property
    def D_expression(self) -> tuple[int, tuple[str, ...]]:
        """Structured exact expression ``1 + sum_i C_phys_i``."""
        return (1, tuple(c.symbol for c in self.jet_constants))

    @property
    def existential_constants_positive(self) -> bool:
        """Positivity is a pinned theorem premise/conclusion, not a numeric evaluation."""
        return True

    @property
    def numeric_D_materialized(self) -> bool:
        return False


@dataclass(frozen=True)
class PhysicalTailPrefactorRow:
    """One derivative row of ``m! * 2^-J * D^m * q^P``."""

    derivative_order: int
    factorial_exact: int
    dyadic_prefactor_exact: Fraction
    D_power: int
    physical_power: Fraction
    target_power: Fraction

    def __post_init__(self) -> None:
        if isinstance(self.derivative_order, bool) or not isinstance(self.derivative_order, int):
            raise TypeError("derivative_order must be an integer")
        if self.derivative_order < 0:
            raise ValueError("derivative_order must be nonnegative")
        if self.factorial_exact != factorial(self.derivative_order):
            raise ValueError("factorial_exact drifted from the pinned chain-rule factor")
        if not isinstance(self.dyadic_prefactor_exact, Fraction):
            raise TypeError("dyadic_prefactor_exact must be an exact Fraction")
        if self.D_power != self.derivative_order:
            raise ValueError("physical-chart D exponent must equal the derivative order")
        if not isinstance(self.physical_power, Fraction) or not isinstance(
            self.target_power, Fraction
        ):
            raise TypeError("tail powers must be exact Fractions")
        if self.physical_power < self.target_power:
            raise ValueError("physical tail power misses the requested target")

    @property
    def target_slack(self) -> Fraction:
        return self.physical_power - self.target_power


@dataclass(frozen=True)
class FinitePhysicalTailPrefactorCertificate:
    """Attach the pinned symbolic physical-chart prefactor to one finite tail request."""

    first_omitted: FiniteFirstOmittedPhysicalTailCertificate
    chart_bound: PhysicalChartFiniteBoundShape
    formal_revision: str = PINNED_FORMAL_REVISION
    composition_theorem: str = PINNED_COMPOSITION_THEOREM
    physical_tail_theorem: str = PINNED_PHYSICAL_TAIL_THEOREM

    def __post_init__(self) -> None:
        if not isinstance(self.first_omitted, FiniteFirstOmittedPhysicalTailCertificate):
            raise TypeError("first_omitted must be a FiniteFirstOmittedPhysicalTailCertificate")
        if not isinstance(self.chart_bound, PhysicalChartFiniteBoundShape):
            raise TypeError("chart_bound must be a PhysicalChartFiniteBoundShape")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.composition_theorem != PINNED_COMPOSITION_THEOREM:
            raise ValueError("composition_theorem does not match physical_composition_bound")
        if self.physical_tail_theorem != PINNED_PHYSICAL_TAIL_THEOREM:
            raise ValueError("physical_tail_theorem does not match exists_physical_uncut_tail")

        M = self.first_omitted.uncut.physical.max_derivative_order
        if self.chart_bound.max_derivative_order != M:
            raise ValueError("physical-chart finite bound uses a different derivative budget")
        if self.chart_bound.h_exact != self.first_omitted.h_exact:
            raise ValueError("physical-chart finite bound uses a different h")

        supported = (
            self.first_omitted.uncut.physical.ordinary.chain.chain.supported_prefix
        )
        x_lo, x_hi, _, _ = supported.exact_support_box
        if (self.chart_bound.x_lower, self.chart_bound.x_upper) != (x_lo, x_hi):
            raise ValueError("physical-chart finite bound is not tied to the hierarchy support box")

        rows = self.rows
        if len(rows) != M + 1:
            raise RuntimeError("physical prefactor ledger has the wrong derivative budget")
        source_rows = self.first_omitted.jet_powers
        if tuple(row.physical_power for row in rows) != tuple(
            row.physical_tail_power for row in source_rows
        ):
            raise ValueError("physical q-powers drifted from the first-omitted ledger")
        if any(
            row.dyadic_prefactor_exact != self.first_omitted.dyadic_prefactor_exact
            for row in rows
        ):
            raise ValueError("physical prefactor lost the common dyadic factor")

    @property
    def rows(self) -> tuple[PhysicalTailPrefactorRow, ...]:
        return tuple(
            PhysicalTailPrefactorRow(
                derivative_order=row.derivative_order,
                factorial_exact=factorial(row.derivative_order),
                dyadic_prefactor_exact=self.first_omitted.dyadic_prefactor_exact,
                D_power=row.derivative_order,
                physical_power=row.physical_tail_power,
                target_power=row.target_power,
            )
            for row in self.first_omitted.jet_powers
        )

    @property
    def truncation_order(self) -> int:
        return self.first_omitted.truncation_order

    @property
    def first_omitted_order(self) -> int:
        return self.first_omitted.first_omitted_order

    @property
    def physical_chart_finite_bound_shape_verified(self) -> bool:
        return True

    @property
    def physical_composition_prefactor_shape_verified(self) -> bool:
        return True

    @property
    def numeric_physical_chart_constant_materialized(self) -> bool:
        return False

    @property
    def actual_pde_residual_verified(self) -> bool:
        return False

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
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


def certify_finite_physical_tail_prefactor(
    h: float,
    max_derivative_order: int,
    target_physical_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FinitePhysicalTailPrefactorCertificate:
    """Build one provider-owned finite tail and attach the pinned chain-rule prefactor shape."""
    first = certify_finite_first_omitted_physical_tail(
        h,
        max_derivative_order=max_derivative_order,
        target_physical_power=target_physical_power,
        provider=provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    supported = first.uncut.physical.ordinary.chain.chain.supported_prefix
    x_lo, x_hi, _, _ = supported.exact_support_box
    chart = PhysicalChartFiniteBoundShape(
        max_derivative_order=max_derivative_order,
        h_exact=first.h_exact,
        x_lower=x_lo,
        x_upper=x_hi,
        jet_constants=tuple(
            PhysicalChartExistentialJetConstant(
                derivative_order=i,
                symbol=f"C_phys_{i}",
            )
            for i in range(max_derivative_order + 1)
        ),
    )
    return FinitePhysicalTailPrefactorCertificate(first_omitted=first, chart_bound=chart)
