from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_residual_rate_target import (
    PINNED_EVENTUALLY_GE_SYMBOL,
    PINNED_GRAPH_LOSS_SYMBOL,
    PINNED_RESIDUAL_RATE_SYMBOL,
    PINNED_TENDSTO_SYMBOL,
    Section9ResidualRateTargetCertificate,
    coordinate_A,
    graph_loss,
    residual_loss,
    residual_rate,
    select_residual_rate_stage,
)


def test_pinned_loss_and_rate_formulas_are_exact():
    h = Fraction(1, 2)
    beta = Fraction(1, 3)

    assert coordinate_A(h) == Fraction(1)
    assert graph_loss(2) == Fraction(8)
    assert residual_loss(h, beta, 2) == Fraction(73, 6)
    assert residual_rate(h, beta, 0, 2) == Fraction(-709, 60)
    assert residual_rate(h, beta, 1, 2) - residual_rate(h, beta, 0, 2) == Fraction(1, 20)


def test_target_selector_returns_the_minimal_exact_stage():
    certificate = select_residual_rate_stage(
        h=Fraction(1, 2),
        beta=Fraction(1, 3),
        derivative_order=2,
        target=5,
    )

    assert certificate.stage == 337
    assert certificate.selected_rate == Fraction(151, 30)
    assert certificate.previous_rate == Fraction(299, 60)
    assert certificate.selected_rate >= certificate.target
    assert certificate.previous_rate < certificate.target
    assert certificate.target_margin == Fraction(1, 30)
    assert certificate.exact_minimal_stage_machine_checked is True


def test_stage_zero_is_selected_when_target_is_already_met():
    certificate = Section9ResidualRateTargetCertificate(
        h=1,
        beta=0,
        derivative_order=0,
        target=-4,
    )

    assert certificate.stage == 0
    assert certificate.selected_rate == Fraction(-19, 5)
    assert certificate.previous_rate is None
    assert certificate.target_margin == Fraction(1, 5)


def test_selector_rejects_approximate_and_invalid_inputs():
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        select_residual_rate_stage(h=0.5, beta=0, derivative_order=0, target=0)
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        select_residual_rate_stage(h=1, beta=Decimal("0.1"), derivative_order=0, target=0)
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        select_residual_rate_stage(h=1, beta=0, derivative_order=0, target=0.1)
    with pytest.raises(ValueError, match="h must be strictly positive"):
        select_residual_rate_stage(h=0, beta=0, derivative_order=0, target=0)
    with pytest.raises(ValueError, match="derivative_order must be nonnegative"):
        select_residual_rate_stage(h=1, beta=0, derivative_order=-1, target=0)
    with pytest.raises(TypeError, match="derivative_order must be a nonnegative integer"):
        select_residual_rate_stage(h=1, beta=0, derivative_order=Fraction(1, 2), target=0)


def test_truth_boundary_and_pinned_dependency_chain_remain_explicit():
    certificate = select_residual_rate_stage(
        h=Fraction(3, 4),
        beta=Fraction(2, 5),
        derivative_order=3,
        target=Fraction(9, 2),
    )

    assert PINNED_GRAPH_LOSS_SYMBOL == "NavierStokes.PhysicalGraphBounds.graphLoss"
    assert PINNED_RESIDUAL_RATE_SYMBOL == "NavierStokes.ActualIterationLedger.residualRate"
    assert PINNED_TENDSTO_SYMBOL.endswith("residualRate_tendsto_atTop")
    assert PINNED_EVENTUALLY_GE_SYMBOL.endswith("eventually_residualRate_ge")
    assert certificate.fixed_derivative_target_only is True
    assert certificate.actual_stage_field_consumed is False
    assert certificate.actual_residual_evaluated is False
    assert certificate.infinite_correction_sequence_certified is False
    assert certificate.eq_9_21_summed_field_certified is False
    assert certificate.paper_exact_velocity_available is False
