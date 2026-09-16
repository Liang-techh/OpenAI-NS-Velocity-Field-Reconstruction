from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_oscillatory_curl_replay import (
    SOURCE_ARCHIVE_PATH,
    SOURCE_ASSOCIATED_CHECKER,
    SOURCE_CHECKER_SHA256,
    SOURCE_COMPONENT,
    SOURCE_COMPONENT_SHA256,
    SOURCE_RESULT_SHA256,
    cartesianized_z_mutation,
    dot,
    normalized_cylindrical_curl_remainder,
    transverse_phase_recovery,
    vector_residual,
)


def test_pinned_public_provenance() -> None:
    assert SOURCE_COMPONENT == "NS-oscillations"
    assert SOURCE_ARCHIVE_PATH == "proof_sources/oscillations/oscillations_body.tex"
    assert SOURCE_ASSOCIATED_CHECKER == "proof_sources/oscillations/exact_checks.py"
    assert SOURCE_COMPONENT_SHA256 == "3fd5c61de39b6f65a581ead5b29c741c2e0f46fb30f62e582dc24d93bbe83615"
    assert SOURCE_CHECKER_SHA256 == "a44ba49990c3b604c883d5b4e1726cb33068867cd85b23f3a2c22847e4ce8184"
    assert SOURCE_RESULT_SHA256 == "30d7a4d934e984a77cfcb04188d72c35391c46e38bcbf25a2024e38bcb783a12"


def test_transverse_phase_polarization_recovers_target_exactly() -> None:
    n = (F(1), F(2), F(-1))
    t = (F(2), F(-1), F(0))
    assert dot(n, t) == 0
    assert transverse_phase_recovery(n, t) == t
    assert vector_residual(transverse_phase_recovery(n, t), t) == (0, 0, 0)


def test_nontransverse_target_fails_closed_as_projection_defect() -> None:
    n = (F(1), F(2), F(-1))
    t = (F(2), F(-1), F(1))
    assert dot(n, t) == -1
    recovered = transverse_phase_recovery(n, t)
    assert recovered != t
    assert vector_residual(recovered, t) == (F(1, 6), F(1, 3), F(-1, 6))


def test_normalized_cylindrical_remainder_fixture_exact() -> None:
    actual = normalized_cylindrical_curl_remainder(
        epsilon=F(1, 8),
        radius=F(3, 2),
        d_z_c_theta=F(2, 5),
        d_theta_c_z=F(7, 9),
        d_z_c_r=F(-3, 8),
        d_r_c_z=F(5, 11),
        c_theta=F(-4, 9),
        d_r_c_theta=F(11, 12),
        d_theta_c_r=F(-2, 7),
    )
    assert actual == (F(-181, 540), F(-19, 44), F(613, 6048))


def test_dropping_cylindrical_product_term_is_detected_exactly() -> None:
    correct_z = normalized_cylindrical_curl_remainder(
        epsilon=F(1, 8),
        radius=F(3, 2),
        d_z_c_theta=F(2, 5),
        d_theta_c_z=F(7, 9),
        d_z_c_r=F(-3, 8),
        d_r_c_z=F(5, 11),
        c_theta=F(-4, 9),
        d_r_c_theta=F(11, 12),
        d_theta_c_r=F(-2, 7),
    )[2]
    wrong_z = cartesianized_z_mutation(
        epsilon=F(1, 8),
        radius=F(3, 2),
        d_r_c_theta=F(11, 12),
        d_theta_c_r=F(-2, 7),
    )
    assert correct_z - wrong_z == F(-1, 27)
    assert correct_z != wrong_z


def test_float_and_zero_radius_inputs_fail_closed() -> None:
    with pytest.raises(TypeError):
        transverse_phase_recovery((1.0, F(2), F(-1)), (F(2), F(-1), F(0)))
    with pytest.raises(ValueError):
        normalized_cylindrical_curl_remainder(
            epsilon=F(1, 8),
            radius=0,
            d_z_c_theta=0,
            d_theta_c_z=0,
            d_z_c_r=0,
            d_r_c_z=0,
            c_theta=0,
            d_r_c_theta=0,
            d_theta_c_r=0,
        )
