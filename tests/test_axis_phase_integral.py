"""Independent checks for the finite validated natural-axis phase integral.

The polynomial-kernel cases below are synthetic exact-rational fixtures.  They
exercise the interval construction and do not certify the paper's schedule
parameters or the eventual amplitude/reconstruction.
"""

from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_phase_integral import (
    NeedSubdivision,
    PhaseIntegrationLimit,
    validated_phase_cell,
    validated_phase_integral,
)


F = Fraction


def _poly_add(left: dict[int, F], right: dict[int, F]) -> dict[int, F]:
    result = dict(left)
    for degree, value in right.items():
        result[degree] = result.get(degree, F(0)) + value
        if result[degree] == 0:
            del result[degree]
    return result


def _poly_scale(poly: dict[int, F], scalar: F) -> dict[int, F]:
    return {degree: scalar * value for degree, value in poly.items() if value}


def _poly_mul(left: dict[int, F], right: dict[int, F]) -> dict[int, F]:
    result: dict[int, F] = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            result[degree] = result.get(degree, F(0)) + left_value * right_value
    return {degree: value for degree, value in result.items() if value}


def _poly_integral(poly: dict[int, F], left: F, right: F) -> F:
    return sum(
        (value * (right ** (degree + 1) - left ** (degree + 1)) / (degree + 1)
         for degree, value in poly.items()),
        F(0),
    )


def test_strict_inputs_zero_and_oriented_phase_integrals() -> None:
    with pytest.raises(TypeError):
        validated_phase_cell(0.0, F(0), F(10), F(0), F(1, 10))
    with pytest.raises(ValueError, match="sigma"):
        validated_phase_cell(F(0), F(0), F(0), F(0), F(1, 10))
    with pytest.raises(ValueError, match="window"):
        validated_phase_integral(
            F(0), F(0), F(10), F(6, 5), absolute_tolerance=F(1, 10**20)
        )
    with pytest.raises(ValueError, match="positive"):
        validated_phase_integral(
            F(0), F(0), F(10), F(1, 10), absolute_tolerance=F(0)
        )

    zero = validated_phase_integral(
        F(0), F(0), F(10), F(0), absolute_tolerance=F(1, 10**20)
    )
    assert zero.integral_estimate == zero.lower == zero.upper == 0
    assert zero.error_bound == 0
    assert zero.cell_count == 0
    assert zero.paper_exact is False

    # Compare the signed run to the independently ordered negative cell.
    tolerance = F(1, 100)
    negative_cell = validated_phase_cell(
        F(0), F(1), F(10), F(-1, 10), F(0), order=16
    )
    negative = validated_phase_integral(
        F(0), F(1), F(10), F(-1, 10),
        absolute_tolerance=tolerance,
        initial_order=16,
        max_order=16,
    )
    assert negative.integral_estimate == -negative_cell.estimate
    assert negative.lower == -negative_cell.upper
    assert negative.upper == -negative_cell.lower


def test_symmetric_odd_cell_is_exactly_zero_and_gate_is_strict() -> None:
    symmetric = validated_phase_cell(
        F(0), F(0), F(10), F(-1, 10), F(1, 10), order=16
    )
    assert symmetric.integral_estimate == 0
    assert symmetric.lower <= 0 <= symmetric.upper
    assert symmetric.lower == -symmetric.upper
    assert 0 <= symmetric.theta < F(1, 2)

    # A small sigma makes the same cell fail the Taylor denominator gate;
    # callers must subdivide instead of accepting an uncertified cell.
    with pytest.raises(NeedSubdivision):
        validated_phase_cell(
            F(0), F(0), F(1), F(-1, 10), F(1, 10), order=16
        )


def test_adaptation_increases_order_and_small_caps_fail_closed() -> None:
    tolerance = F(1, 10**24)
    result = validated_phase_integral(
        F(0), F(0), F(10), F(1, 10),
        absolute_tolerance=tolerance,
        initial_order=2,
        max_order=256,
    )
    assert result.error_bound <= tolerance
    assert result.max_order_used > 2
    assert result.cell_count > 0

    with pytest.raises(PhaseIntegrationLimit):
        validated_phase_integral(
            F(0), F(0), F(10), F(1, 10),
            absolute_tolerance=tolerance,
            initial_order=2,
            max_order=2,
        )


def test_adaptive_subdivision_and_geometry_caps_fail_closed() -> None:
    # This synthetic sigma=1 interval fails the initial denominator gate and
    # then succeeds after genuine bisection into multiple accepted cells.
    kwargs = dict(
        h=F(0),
        j=F(0),
        sigma=F(1),
        eta=F(1, 2),
        absolute_tolerance=F(1, 10**10),
        initial_order=16,
        max_order=128,
    )
    result = validated_phase_integral(**kwargs, max_cells=4096, max_depth=128)
    assert result.cell_count > 1
    assert result.error_bound <= kwargs["absolute_tolerance"]

    with pytest.raises(PhaseIntegrationLimit, match="cell"):
        validated_phase_integral(**kwargs, max_cells=1, max_depth=128)
    with pytest.raises(PhaseIntegrationLimit, match="depth"):
        validated_phase_integral(**kwargs, max_cells=4096, max_depth=1)


def test_independent_global_geometric_oracle_overlaps_validated_interval() -> None:
    # Synthetic exact-rational kernel: H(x) = 9*x/2 - 4*x^3, L(x) = 1,
    # sigma = 10.  Expand 1/(sigma^2 + H^2) globally in H^2/sigma^2;
    # this is independent of the cell-center Taylor recurrence.
    eta = F(1, 10)
    sigma = F(10)
    H = {1: F(9, 2), 3: F(-4)}
    H_squared_over_sigma_squared = _poly_scale(
        _poly_mul(H, H), F(1, 100)
    )
    terms = {}
    term = {0: F(1)}
    truncation_order = 10
    for k in range(truncation_order + 1):
        terms = _poly_add(terms, _poly_scale(term, F((-1) ** k)))
        term = _poly_mul(term, H_squared_over_sigma_squared)
    partial_integrand = _poly_scale(_poly_mul(H, terms), F(-1, 100))
    partial = _poly_integral(partial_integrand, F(0), eta)

    h_bound = sum(
        (abs(value) * eta**degree for degree, value in H.items()),
        F(0),
    )
    q = h_bound * h_bound / sigma**2
    assert q < 1
    residual = eta * h_bound / sigma**2 * q ** (truncation_order + 1) / (1 - q)
    oracle_lower = partial - residual
    oracle_upper = partial + residual

    tolerance = F(1, 10**36)
    result = validated_phase_integral(
        F(0), F(0), sigma, eta,
        absolute_tolerance=tolerance,
        initial_order=16,
        max_order=512,
    )
    assert result.lower <= oracle_upper
    assert oracle_lower <= result.upper
    assert result.upper - result.lower < oracle_upper - oracle_lower
    assert result.error_bound <= tolerance
