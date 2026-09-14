from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_finite_jet_residual_target import (
    PINNED_EVENTUALLY_GE_SYMBOL,
    PINNED_LEDGER_FILE,
    PINNED_RESIDUAL_RATE_SYMBOL,
    Section9FiniteJetResidualTargetCertificate,
    select_finite_jet_residual_stage,
)


def test_one_common_stage_covers_every_finite_derivative_order_exactly():
    certificate = select_finite_jet_residual_stage(
        h=Fraction(1, 2),
        beta=Fraction(1, 3),
        max_derivative=3,
        target=5,
    )

    assert certificate.derivative_orders == (0, 1, 2, 3)
    assert certificate.per_order_minimal_stages == (163, 230, 337, 483)
    assert certificate.stage == 483
    assert certificate.active_orders == (3,)
    assert certificate.selected_rates == (
        Fraction(21),
        Fraction(53, 3),
        Fraction(37, 3),
        Fraction(5),
    )
    assert certificate.previous_rates == (
        Fraction(419, 20),
        Fraction(1057, 60),
        Fraction(737, 60),
        Fraction(99, 20),
    )
    assert all(rate >= certificate.target for rate in certificate.selected_rates)
    assert any(rate < certificate.target for rate in certificate.previous_rates)
    assert certificate.joint_target_margin == 0
    assert certificate.one_common_prefix_machine_checked is True


def test_selector_does_not_assume_derivative_order_monotonicity():
    # A sufficiently negative exact beta can make a lower derivative order the
    # active constraint.  The implementation must therefore inspect every
    # finite row instead of silently selecting only m=M.
    certificate = select_finite_jet_residual_stage(
        h=1,
        beta=-10,
        max_derivative=3,
        target=0,
    )

    assert certificate.per_order_minimal_stages == (38, 0, 0, 0)
    assert certificate.stage == 38
    assert certificate.active_orders == (0,)
    assert certificate.selected_rates[0] == 0
    assert certificate.previous_rates[0] == Fraction(-1, 10)


def test_stage_zero_is_valid_only_when_all_requested_orders_already_pass():
    certificate = Section9FiniteJetResidualTargetCertificate(
        h=1,
        beta=0,
        max_derivative=2,
        target=-20,
    )

    assert certificate.stage == 0
    assert certificate.previous_rates == ()
    assert all(rate >= certificate.target for rate in certificate.selected_rates)


def test_joint_selector_rejects_approximate_or_invalid_theorem_inputs():
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        select_finite_jet_residual_stage(h=0.5, beta=0, max_derivative=2, target=0)
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        select_finite_jet_residual_stage(h=1, beta=Decimal("0.1"), max_derivative=2, target=0)
    with pytest.raises(TypeError, match="float/Decimal approximations are rejected"):
        select_finite_jet_residual_stage(h=1, beta=0, max_derivative=2, target=0.1)
    with pytest.raises(ValueError, match="max_derivative must be nonnegative"):
        select_finite_jet_residual_stage(h=1, beta=0, max_derivative=-1, target=0)
    with pytest.raises(TypeError, match="max_derivative must be a nonnegative integer"):
        select_finite_jet_residual_stage(
            h=1,
            beta=0,
            max_derivative=Fraction(3, 2),
            target=0,
        )


def test_truth_boundary_and_pinned_dependencies_are_explicit():
    certificate = select_finite_jet_residual_stage(
        h=Fraction(1, 4),
        beta=Fraction(2, 5),
        max_derivative=4,
        target=Fraction(7, 3),
    )

    assert PINNED_LEDGER_FILE == "NavierStokes/ActualIterationLedger.lean"
    assert PINNED_RESIDUAL_RATE_SYMBOL.endswith("ActualIterationLedger.residualRate")
    assert PINNED_EVENTUALLY_GE_SYMBOL.endswith("eventually_residualRate_ge")
    assert certificate.finite_derivative_budget_only is True
    assert certificate.actual_stage_field_consumed is False
    assert certificate.actual_residual_evaluated is False
    assert certificate.all_derivative_orders_certified is False
    assert certificate.infinite_correction_sequence_certified is False
    assert certificate.eq_9_21_summed_field_certified is False
    assert certificate.paper_exact_velocity_available is False
