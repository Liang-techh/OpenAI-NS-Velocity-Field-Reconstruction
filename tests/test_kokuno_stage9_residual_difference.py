from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_stage9_residual_difference import (
    componentwise_a39_majorant_m0,
    expanded_difference_zero_order,
    replay_zero_order_difference,
    verify_componentwise_majorant,
)


FIXTURE = dict(
    v=(F(2, 3), F(-1, 4), F(5, 6)),
    e=(F(1, 5), F(-2, 7), F(3, 8)),
    dt_v=(F(1, 9), F(-2, 11), F(3, 13)),
    dt_e=(F(-4, 15), F(5, 17), F(-6, 19)),
    grad_v=(
        (F(1, 2), F(-1, 3), F(2, 5)),
        (F(3, 7), F(1, 4), F(-2, 9)),
        (F(-1, 6), F(5, 12), F(7, 15)),
    ),
    grad_e=(
        (F(-2, 11), F(3, 10), F(1, 7)),
        (F(4, 13), F(-5, 14), F(2, 9)),
        (F(1, 8), F(-3, 16), F(5, 18)),
    ),
    lap_v=(F(7, 20), F(-3, 10), F(11, 30)),
    lap_e=(F(-5, 21), F(4, 25), F(-7, 27)),
    grad_p=(F(2, 9), F(-1, 8), F(5, 22)),
    grad_pi=(F(-3, 14), F(7, 24), F(-2, 15)),
)


def test_direct_residual_difference_matches_complete_expansion_exactly():
    replay = replay_zero_order_difference(**FIXTURE)
    assert replay.exact
    assert replay.residual == (F(0), F(0), F(0))
    assert replay.direct == replay.expanded


def test_fixture_has_pinned_exact_nonzero_difference():
    replay = replay_zero_order_difference(**FIXTURE)
    assert replay.direct == (
        F(-10, 231),
        F(1218167, 1124550),
        F(433429, 1149120),
    )


def test_omitting_quadratic_error_transport_is_detected():
    full = expanded_difference_zero_order(
        v=FIXTURE["v"],
        e=FIXTURE["e"],
        dt_e=FIXTURE["dt_e"],
        grad_v=FIXTURE["grad_v"],
        grad_e=FIXTURE["grad_e"],
        lap_e=FIXTURE["lap_e"],
        grad_pi=FIXTURE["grad_pi"],
    )
    # Delete e_k * partial_k e_i only, retaining all other A38 terms.
    wrong = []
    for i in range(3):
        quadratic = sum(FIXTURE["e"][k] * FIXTURE["grad_e"][i][k] for k in range(3))
        wrong.append(full[i] - quadratic)
    defect = tuple(full[i] - wrong[i] for i in range(3))
    assert defect == (F(-211, 3080), F(9437, 38220), F(307, 1680))
    assert defect != (F(0), F(0), F(0))


def test_pressure_increment_sign_mutation_fails_closed():
    replay = replay_zero_order_difference(**FIXTURE)
    wrong = list(replay.expanded)
    # Flip +partial_1 pi to -partial_1 pi in component 1.
    wrong[0] -= 2 * FIXTURE["grad_pi"][0]
    assert replay.direct[0] - wrong[0] == 2 * FIXTURE["grad_pi"][0]
    assert replay.direct[0] != wrong[0]


def test_componentwise_a39_m0_majorant_is_exact_rational():
    v2 = F(5, 3)
    e2 = F(7, 5)
    bound = componentwise_a39_majorant_m0(v2_bound=v2, e2_bound=e2)
    assert bound == F(672, 25)
    assert verify_componentwise_majorant(
        difference=replay_zero_order_difference(**FIXTURE).direct,
        majorant=bound,
    )


def test_zero_error_gives_zero_difference_and_zero_majorant():
    zero = (F(0), F(0), F(0))
    zero_grad = (zero, zero, zero)
    data = dict(FIXTURE)
    data.update(e=zero, dt_e=zero, grad_e=zero_grad, lap_e=zero, grad_pi=zero)
    replay = replay_zero_order_difference(**data)
    assert replay.direct == zero
    assert componentwise_a39_majorant_m0(v2_bound=F(17, 4), e2_bound=F(0)) == 0


def test_approximate_and_negative_bound_inputs_fail_closed():
    with pytest.raises(TypeError):
        componentwise_a39_majorant_m0(v2_bound=0.5, e2_bound=F(1))  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        componentwise_a39_majorant_m0(v2_bound=F(1), e2_bound=F(-1, 10))
    bad = dict(FIXTURE)
    bad["e"] = (0.2, F(-2, 7), F(3, 8))
    with pytest.raises(TypeError):
        replay_zero_order_difference(**bad)
