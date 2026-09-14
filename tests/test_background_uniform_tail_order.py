from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_uniform_tail_order import (
    PINNED_DIAGONAL_THEOREM,
    PINNED_FORMAL_REVISION,
    PINNED_SLOWBOREL_THEOREM,
    PinnedUniformTailOrderCertificate,
    certify_pinned_uniform_tail_order,
)


def test_selects_exact_minimal_pinned_tail_order():
    cert = certify_pinned_uniform_tail_order(
        0.01,
        Fraction(1, 2),
        Fraction(3, 1),
        minimum_order=7,
    )

    assert cert.tail_power_exact >= Fraction(3, 1)
    assert cert.gain_next_exact >= 2 * (Fraction(3, 1) + Fraction(1, 2))
    assert cert.target_margin_exact >= 0
    assert cert.dyadic_prefactor_exact == Fraction(1, 2) ** cert.order

    # The target, not minimum_order=7, determines this large J.  The immediately
    # preceding stage must therefore fail the exact requested power.
    assert cert.order > cert.minimum_order
    assert cert.predecessor_fails_target_when_not_minimum_limited
    previous_power = cert.h_exact * cert.order - cert.derivative_loss
    assert previous_power < cert.target_power


def test_minimum_order_can_dominate_target_without_changing_theorem_inequality():
    cert = certify_pinned_uniform_tail_order(
        0.01,
        0,
        Fraction(1, 100),
        minimum_order=25,
    )
    assert cert.order == 25
    assert cert.tail_power_exact >= cert.target_power
    assert not cert.predecessor_fails_target_when_not_minimum_limited


@pytest.mark.parametrize("target", [1, Fraction(7, 3), 25])
def test_exact_selector_matches_bruteforce_minimum(target):
    loss = Fraction(2, 5)
    minimum = 4
    cert = certify_pinned_uniform_tail_order(0.01, loss, target, minimum_order=minimum)

    admissible = [
        J
        for J in range(minimum, cert.order + 1)
        if cert.h_exact * (J + 1) - loss >= Fraction(target)
    ]
    assert admissible
    assert cert.order == min(admissible)


def test_certificate_replays_pinned_source_identity_and_stays_fail_closed():
    cert = certify_pinned_uniform_tail_order(0.01, 2, 5)
    assert cert.formal_revision == PINNED_FORMAL_REVISION
    assert cert.diagonal_theorem == PINNED_DIAGONAL_THEOREM
    assert cert.slowborel_theorem == PINNED_SLOWBOREL_THEOREM
    assert cert.pinned_scalar_tail_arithmetic

    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_float_target_or_loss_is_rejected_instead_of_rounding_a_theorem_gate():
    with pytest.raises(TypeError):
        certify_pinned_uniform_tail_order(0.01, 0.5, 3)
    with pytest.raises(TypeError):
        certify_pinned_uniform_tail_order(0.01, 0, 3.0)


def test_invalid_exact_parameters_fail_closed():
    with pytest.raises(ValueError):
        certify_pinned_uniform_tail_order(0.01, -1, 3)
    with pytest.raises(ValueError):
        certify_pinned_uniform_tail_order(0.01, 0, 0)
    with pytest.raises(ValueError):
        certify_pinned_uniform_tail_order(0.01, 0, 1, minimum_order=-1)


def test_constructor_rejects_nonminimal_or_cross_revision_witnesses():
    good = certify_pinned_uniform_tail_order(0.01, 1, 4, minimum_order=2)

    with pytest.raises(ValueError):
        PinnedUniformTailOrderCertificate(
            h=good.h,
            derivative_loss=good.derivative_loss,
            target_power=good.target_power,
            minimum_order=good.minimum_order,
            order=good.order + 1,
        )

    with pytest.raises(ValueError):
        PinnedUniformTailOrderCertificate(
            h=good.h,
            derivative_loss=good.derivative_loss,
            target_power=good.target_power,
            minimum_order=good.minimum_order,
            order=good.order,
            formal_revision="wrong-revision",
        )
