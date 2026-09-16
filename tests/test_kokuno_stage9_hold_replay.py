from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_stage9_hold_replay import (
    SOURCE_BUNDLE_SHA256,
    SOURCE_RECORD_COMMIT,
    STAGE9_BODY_SHA256,
    STAGE9_CHECKER_PATH,
    concrete_endpoint_ratios,
    replay_stage9_l_minus_one_hold,
)


def test_stage9_hold_scaling_replay_is_exact() -> None:
    result = replay_stage9_l_minus_one_hold()

    assert result.tolerance == Fraction(0)
    assert result.l_value == Fraction(-1)
    assert result.xe2_log_rate == Fraction(-2)
    assert result.e_log_rate == Fraction(-3, 2)
    assert result.xe2_endpoint_h_power == Fraction(8)
    assert result.e_endpoint_h_power == Fraction(6)
    assert result.rates_exact
    assert result.endpoint_powers_exact
    assert result.passed


def test_stage9_hold_scaling_matches_exact_rational_fixture() -> None:
    h = Fraction(1, 16)
    xe2_ratio, e_ratio = concrete_endpoint_ratios(h)

    assert xe2_ratio == Fraction(1, 16**8)
    assert e_ratio == Fraction(1, 16**6)


def test_stage9_hold_scaling_rejects_length_drift() -> None:
    drifted = replay_stage9_l_minus_one_hold(
        hold_length_log_inverse_h_coefficient=Fraction(4) + Fraction(1, 2**40)
    )

    assert drifted.tolerance == 0
    assert drifted.rates_exact
    assert not drifted.endpoint_powers_exact
    assert not drifted.passed


def test_stage9_hold_scaling_rejects_out_of_domain_h() -> None:
    with pytest.raises(ValueError):
        concrete_endpoint_ratios(Fraction(1))


def test_stage9_source_provenance_is_pinned() -> None:
    assert SOURCE_RECORD_COMMIT == "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
    assert SOURCE_BUNDLE_SHA256 == "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
    assert STAGE9_BODY_SHA256 == "bce636fb3ed7348ed9c73912183f81a75557cc975418ce9e162c4be994f82cb0"
    assert STAGE9_CHECKER_PATH == "proof_sources/stage9/exact_checks.py"
