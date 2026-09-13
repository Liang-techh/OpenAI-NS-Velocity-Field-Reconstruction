from decimal import Decimal, localcontext
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_mixed_scale import (
    MixedScaleCoefficient,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_profile_budget import (
    actual_schedule_profile_budget,
)
from openai_ns_reconstruction.natural_axis_wide import wide_natural_profile_prefix
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_axis_pressure import axis_pressure


def _data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def solver():
    return formal_axis_coefficient_solver(
        actual_schedule_wide_first_picard_state(_data(), 0.05)
    )


def _channel(coefficient: MixedScaleCoefficient, key: tuple[int, int]) -> Decimal:
    return coefficient.channels.get(key, Decimal(0))


def test_pressure_primitive_and_sparse_channels_are_independent(solver) -> None:
    eta = 0.17
    rows = solver.profile_jet_prefix(2, 0, eta)
    pressure0, pressure1, pressure2 = (row[2] for row in rows)

    # P = primitive(a^2 phi^2): P0=0, P1=a^2 phi0^2, and
    # P2=a^2*(2 phi0 phi1)/2.  Build the last relation by hand.
    assert pressure0 == MixedScaleCoefficient.zero()
    assert pressure1 == MixedScaleCoefficient.channel(2, 0, Decimal(1))
    phi1 = solver.jet_pair(1, 0, eta)[0]
    expected_pressure2 = MixedScaleCoefficient(
        {(q + 2, p): value for (q, p), value in phi1.items()}
    )
    assert pressure2 == expected_pressure2

    odd = MixedScaleCoefficient.channel(1, 0, Decimal(1))
    assert odd.support == ((1, 0),)
    with pytest.raises(ValueError, match="nonnegative even|nonnegative"):
        MixedScaleCoefficient.channel(-1, 0, Decimal(1))


def test_wide_profile_at_axis_keeps_direct_fields_and_false_flags(solver) -> None:
    eta = 0.17
    result = wide_natural_profile_prefix(solver, 2, Decimal(0), eta)
    reference = solver.x1.reference.reference
    with localcontext() as context:
        context.prec = 96
        axis_u = +(Decimal(4) * Decimal.from_float(eta) + Decimal.from_float(reference.j))
    expected_pressure = Decimal.from_float(axis_pressure(_data(), eta))

    assert result.E == MixedScaleCoefficient.zero()
    assert result.F == MixedScaleCoefficient.channel(1, 0, Decimal(1))
    assert result.U == MixedScaleCoefficient.channel(0, 0, axis_u)
    assert result.average_U == result.U
    assert result.dU_deta == MixedScaleCoefficient.channel(0, 0, Decimal(4))
    assert result.d_average_U_deta == MixedScaleCoefficient.channel(0, 0, Decimal(4))
    assert result.V0 == MixedScaleCoefficient.zero()
    assert result.Pi == MixedScaleCoefficient.channel(0, 0, expected_pressure)
    assert result.global_axis_norm_certified is False
    assert result.truncation_certified is False
    assert result.paper_exact is False


def test_wide_profile_at_nonzero_x_preserves_odd_logs_and_pressure_shift(solver) -> None:
    eta = 0.17
    with localcontext() as context:
        context.prec = 96
        X = +(Decimal(1) / solver.Lambda)
        sqrt_factor = +(Decimal(2) * X).sqrt()
    result = wide_natural_profile_prefix(solver, 2, X, eta)

    assert Fraction(result.Y) == Fraction(solver.Lambda) * Fraction(X)
    assert result.F.support
    assert all(q % 2 == 1 for q, _ in result.F.support)
    assert _channel(result.F, (1, 0)) != 0
    assert result.E == result.F.scale(sqrt_factor)

    F_logs = result.terms_log("F")
    assert F_logs
    assert all(term.sign != 0 for term in F_logs.values())
    assert all(term.log_scale is not None and term.log_factor is not None for term in F_logs.values())

    expected_axis_pressure = Decimal.from_float(axis_pressure(_data(), eta))
    assert _channel(result.Pi, (0, 0)) == expected_axis_pressure
    assert _channel(result.Pi, (2, 1)) != 0
    assert result.terms_log("Pi")

    budget = actual_schedule_profile_budget(solver)
    assert Fraction(budget.angular_ratio_norm_upper) >= (
        Fraction(64)
        * Fraction(budget.amplitude_norm_upper)
        * Fraction(budget.profile_norm_upper)
    )
    assert Fraction(budget.pressure_norm_upper) >= (
        Fraction(80 * 64) * Fraction(budget.angular_ratio_norm_upper) ** 2
    )
    assert budget.conditional_on_axis_space_identification is True
    assert budget.global_axis_norm_certified is False
    assert budget.includes_coefficient_roundoff is False
    assert budget.paper_exact is False


def test_n0_v0_uses_independent_ordinary_polynomial(solver) -> None:
    eta = 0.17
    with localcontext() as context:
        context.prec = 96
        X = +(Decimal(1) / solver.Lambda)
        eta_decimal = Decimal.from_float(eta)
        h = Decimal.from_float(solver.data.h)
        D = Decimal.from_float(solver.data.D)
        L = +(Decimal(1) - Decimal(2) * h * eta_decimal * eta_decimal)
        U = +(Decimal(4) * eta_decimal + Decimal.from_float(solver.x1.reference.reference.j))
        expected = +(
            (X / L)
            * (
                Decimal(2) * eta_decimal * U
                - Decimal(2) * D * eta_decimal * U
                - (Decimal(1) - eta_decimal * eta_decimal) * Decimal(4)
            )
        )
    result = wide_natural_profile_prefix(solver, 0, X, eta)
    assert result.U == MixedScaleCoefficient.channel(0, 0, U)
    assert result.average_U == result.U
    assert result.d_average_U_deta == MixedScaleCoefficient.channel(0, 0, Decimal(4))
    actual = _channel(result.V0, (0, 0))
    scale = max(abs(actual), abs(expected), Decimal("1e-1000"))
    assert abs(actual - expected) <= Decimal("1e-90") * scale


def test_conditional_truncation_bounds_scale_all_fields_upward(solver) -> None:
    eta = 0.17
    with localcontext() as context:
        context.prec = 96
        X = +(Decimal(1) / solver.Lambda)
    profile = wide_natural_profile_prefix(solver, 2, X, eta)
    budget = actual_schedule_profile_budget(solver)
    bounds = profile.conditional_truncation_bounds(budget)
    Y = profile.Y

    f_tail = budget.conditional_angular_ratio_tail_bound(Y, profile.max_n)
    u_tail = budget.conditional_tail_bound(Y, profile.max_n)
    eta_tail = budget.conditional_tail_bound(Y, profile.max_n, eta_order=1)
    average_tail = budget.conditional_tail_bound(Y, profile.max_n, average=True)
    average_eta_tail = budget.conditional_tail_bound(
        Y,
        profile.max_n,
        eta_order=1,
        average=True,
    )
    pressure_tail = budget.conditional_pressure_tail_bound(Y, profile.max_n)
    assert bounds["F"] == f_tail.upper_bound
    for name, tail in (
        ("U", u_tail),
        ("dU_deta", eta_tail),
        ("average_U", average_tail),
        ("d_average_U_deta", average_eta_tail),
        ("Pi", pressure_tail),
    ):
        assert Fraction(bounds[name]) >= (
            Fraction(tail.upper_bound) / Fraction(profile.Lambda)
        )
    assert Fraction(bounds["E"]) ** 2 >= (
        Fraction(2) * Fraction(profile.X) * Fraction(bounds["F"]) ** 2
    )
    eta_fraction = Fraction(Decimal.from_float(eta))
    D_fraction = Fraction(Decimal.from_float(solver.data.D))
    h_fraction = Fraction(Decimal.from_float(solver.data.h))
    L_fraction = Fraction(1) - Fraction(2) * h_fraction * eta_fraction**2
    u_fraction = Fraction(u_tail.upper_bound) / Fraction(profile.Lambda)
    average_fraction = Fraction(average_tail.upper_bound) / Fraction(profile.Lambda)
    average_eta_fraction = Fraction(average_eta_tail.upper_bound) / Fraction(profile.Lambda)
    v0_triangle = abs(Fraction(profile.X) / L_fraction) * (
        abs(Fraction(2) * eta_fraction) * u_fraction
        + abs(Fraction(2) * D_fraction * eta_fraction) * average_fraction
        + abs(Fraction(1) - eta_fraction**2) * average_eta_fraction
    )
    assert Fraction(bounds["V0"]) >= v0_triangle
    assert bounds.conditional_on_axis_space_identification is True
    assert bounds.includes_coefficient_roundoff is False
    assert bounds.includes_axis_pressure_quadrature_roundoff is False
    assert bounds.paper_exact is False

    zero_profile = wide_natural_profile_prefix(solver, 0, Decimal(0), eta)
    zero_bounds = zero_profile.conditional_truncation_bounds(budget)
    assert all(value == Decimal(0) for value in zero_bounds.values())

    other_solver = formal_axis_coefficient_solver(solver.x1)
    other_budget = actual_schedule_profile_budget(other_solver)
    with pytest.raises(ValueError, match="exact solver"):
        profile.conditional_truncation_bounds(other_budget)
