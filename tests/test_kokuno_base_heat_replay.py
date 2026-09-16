from fractions import Fraction

from openai_ns_reconstruction.kokuno_base_heat_replay import (
    BASE_HEAT_CHECKER_SHA256,
    BASE_HEAT_DERIVATION_SHA256,
    BASE_HEAT_RECORDED_RESULT_SHA256,
    SOURCE_BUNDLE_SHA256,
    SOURCE_RECORD_COMMIT,
    replay_base_heat_ode_polynomial,
    same_polynomial,
)


def test_base_heat_ode_polynomial_replay_is_exact() -> None:
    result = replay_base_heat_ode_polynomial()

    assert result.tolerance == Fraction(0)
    assert result.ode_reduction_exact
    assert result.total_derivative_factor_exact
    assert result.radial_heat_coefficient_exact
    assert result.passed

    expected = {
        (0, 0, 0): Fraction(1),
        (1, 0, 0): Fraction(1),
        (0, 0, 1): Fraction(-1),
        (0, 1, 2): Fraction(-1),
    }
    assert dict(result.reduced_polynomial) == expected
    assert dict(result.ode_common_factor_polynomial) == expected
    assert dict(result.total_derivative_factor_polynomial) == expected


def test_base_heat_replay_detects_exact_coefficient_drift() -> None:
    result = replay_base_heat_ode_polynomial()
    drifted = dict(result.reduced_polynomial)
    drifted[(0, 0, 0)] += Fraction(1, 2**40)

    assert not same_polynomial(result.ode_common_factor_polynomial, drifted)


def test_base_heat_source_identities_are_pinned() -> None:
    assert SOURCE_RECORD_COMMIT == "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
    assert SOURCE_BUNDLE_SHA256 == "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
    assert BASE_HEAT_DERIVATION_SHA256 == "98d101542a7d8d5c1475764ec99178a663fe363e416d572eeafa048e724a193d"
    assert BASE_HEAT_CHECKER_SHA256 == "da6fe495af223772824b2f658c83c4aec08d7aa41882818589eaf14c6a3425dc"
    assert BASE_HEAT_RECORDED_RESULT_SHA256 == "5b74d2d2598841ce2512790dd70efe35640b69a4ea5a8516a4dbbd482f0bd843"
