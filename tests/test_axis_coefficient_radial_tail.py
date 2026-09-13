from decimal import Decimal, localcontext
from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_profile_budget import (
    actual_schedule_profile_budget,
)
from openai_ns_reconstruction.axis_coefficient_radial_tail import (
    axis_coefficient_radial_tail_bound,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _finite_omitted_subtotal(
    norm_upper: Decimal,
    epsilon: Decimal,
    Y: Decimal,
    max_n: int,
    eta_order: int,
    radial_order: int,
    average: bool,
    stop: int,
) -> Fraction:
    """Directly sum finitely many omitted weighted terms.

    This subtotal is a lower-bound sanity check for the returned majorant; it
    says nothing about the uncomputed infinite tail beyond ``stop``.
    """

    total = Fraction(0)
    Y_fraction = abs(Fraction(Y))
    first = max(max_n + 1, radial_order)
    for n in range(first, stop + 1):
        weight = (
            Fraction(1, 20**n)
            * Fraction(epsilon) ** (-eta_order)
            * Fraction(math.factorial(eta_order) * math.comb(n + eta_order, eta_order),
                       (n + 1) ** 2 * (eta_order + 1) ** 2)
        )
        radial_factor = (
            Y_fraction ** (n - radial_order)
            * Fraction(math.factorial(n), math.factorial(n - radial_order))
        )
        if average:
            radial_factor /= n + 1
        total += Fraction(norm_upper) * weight * radial_factor
    return total


def test_returned_bound_dominates_independent_finite_subtotal() -> None:
    norm_upper = Decimal("3")
    epsilon = Decimal("0.2")
    Y = Decimal("0.7")
    result = axis_coefficient_radial_tail_bound(
        norm_upper,
        epsilon,
        Y,
        max_n=2,
        eta_order=2,
        radial_order=1,
        average=True,
    )
    subtotal = _finite_omitted_subtotal(
        norm_upper,
        epsilon,
        Y,
        max_n=2,
        eta_order=2,
        radial_order=1,
        average=True,
        stop=30,
    )
    assert Fraction(result.upper_bound) >= subtotal
    assert result.conditional_on_global_norm is True
    assert result.includes_coefficient_roundoff is False
    assert result.paper_exact is False

    precise_Y = Decimal("0.700000000000000000000000000000123456789")
    with localcontext() as ctx:
        ctx.prec = 8
        low_precision_context = axis_coefficient_radial_tail_bound(
            norm_upper,
            epsilon,
            precise_Y,
            max_n=2,
            eta_order=2,
            radial_order=1,
            average=True,
        )
    with localcontext() as ctx:
        ctx.prec = 120
        high_precision_context = axis_coefficient_radial_tail_bound(
            norm_upper,
            epsilon,
            precise_Y,
            max_n=2,
            eta_order=2,
            radial_order=1,
            average=True,
        )
    assert low_precision_context.upper_bound == high_precision_context.upper_bound


def test_y_zero_selects_exact_first_omitted_derivative_row() -> None:
    norm_upper = Decimal("7")
    epsilon = Decimal("0.25")
    eta_order = 2
    radial_order = 5
    max_n = 4
    result = axis_coefficient_radial_tail_bound(
        norm_upper,
        epsilon,
        Decimal(0),
        max_n=max_n,
        eta_order=eta_order,
        radial_order=radial_order,
        average=True,
    )
    n = radial_order
    expected = (
        Fraction(norm_upper)
        * Fraction(1, 20**n)
        * Fraction(epsilon) ** (-eta_order)
        * Fraction(math.factorial(eta_order) * math.comb(n + eta_order, eta_order),
                   (n + 1) ** 2 * (eta_order + 1) ** 2)
        * Fraction(math.factorial(n), math.factorial(n - radial_order))
        / (n + 1)
    )
    assert Fraction(result.upper_bound) >= expected

    no_omitted_derivative = axis_coefficient_radial_tail_bound(
        norm_upper,
        epsilon,
        Decimal(0),
        max_n=5,
        eta_order=eta_order,
        radial_order=2,
    )
    assert no_omitted_derivative.upper_bound == Decimal(0)


def test_radial_tail_domain_guards() -> None:
    with pytest.raises(ValueError):
        axis_coefficient_radial_tail_bound(Decimal(1), Decimal("0.2"), Decimal(20), 2)
    with pytest.raises(ValueError):
        axis_coefficient_radial_tail_bound(Decimal(1), Decimal(0), Decimal("0.1"), 2)
    with pytest.raises(ValueError):
        axis_coefficient_radial_tail_bound(Decimal(1), Decimal("0.2"), Decimal("0.1"), -1)


def test_actual_profile_budget_binds_solver_and_conditional_tail() -> None:
    data = TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )
    x1 = actual_schedule_wide_first_picard_state(data, 0.05)
    solver = formal_axis_coefficient_solver(x1)
    budget = actual_schedule_profile_budget(solver)

    assert budget.Lambda == solver.Lambda
    assert budget.epsilon == Decimal.from_float(solver.epsilon)
    assert Fraction(budget.profile_norm_upper) >= (
        Fraction(budget.reference_norm_upper)
        + Fraction(budget.certificate.one_step_radius_upper)
    )
    assert budget.conditional_on_axis_space_identification is True
    assert budget.global_axis_norm_certified is False
    assert budget.includes_coefficient_roundoff is False
    assert budget.paper_exact is False

    zero_tail = budget.conditional_tail_bound(Decimal(0), max_n=0)
    assert zero_tail.upper_bound == Decimal(0)
    assert zero_tail.conditional_on_global_norm is True
    assert zero_tail.includes_coefficient_roundoff is False
