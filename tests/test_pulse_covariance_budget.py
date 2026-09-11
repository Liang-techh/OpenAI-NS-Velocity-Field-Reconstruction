import math

import numpy as np
import pytest

from openai_ns_reconstruction.pulse_covariance_budget import (
    SignedPulseCovarianceErrorBudget,
    sqrt_scale_error_envelope,
    uniform_two_sign_error_bound,
)
from openai_ns_reconstruction.stress_cone import normal_magnitude


def _weighted_average(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    # Independent oracle used only in tests; production does not evaluate the
    # pulse integral or infer an analytic bound from samples.
    return np.sum(values * weights[:, None], axis=0) / float(np.sum(weights))


def test_reference_direction_matches_signed_stress_cone_convention():
    c_star = -0.65
    u_star = 1.7
    a = normal_magnitude(c_star, u_star)

    sigma_plus = SignedPulseCovarianceErrorBudget(
        sigma=1,
        c_star=c_star,
        u_star=u_star,
        ratio_error_bound=0.0,
        normalized_first_moment_bound=0.0,
    )
    sigma_minus = SignedPulseCovarianceErrorBudget(
        sigma=-1,
        c_star=c_star,
        u_star=u_star,
        ratio_error_bound=0.0,
        normalized_first_moment_bound=0.0,
    )

    # Eq. (7.28): H_sigma/h_sigma = -A_c N - sigma u_* K + e_sigma.
    assert sigma_plus.reference_direction == pytest.approx((-a, -u_star))
    assert sigma_minus.reference_direction == pytest.approx((-a, u_star))
    assert sigma_plus.ideal_direction(u_star) == pytest.approx((-a, -u_star))
    assert sigma_minus.ideal_direction(-u_star) == pytest.approx((-a, u_star))


def test_weighted_pulse_average_obeys_independent_triangle_budget():
    sigma = 1
    c_star = -0.72
    u_star = 1.4
    length = 80.0
    midpoint = length / 2.0

    # Positive synthetic weights stand in only for the algebra of the Eq. (7.27)
    # normalized average.  They are not paper pulse data.
    v = np.array([8.0, 22.0, 35.0, 40.0, 47.0, 61.0, 76.0])
    weights = np.array([0.2, 0.8, 2.5, 4.0, 2.1, 0.7, 0.15])
    s = sigma * (u_star / 2.0 + u_star * v / length)

    ideal = np.column_stack(
        (
            c_star * np.sqrt(1.0 + np.square(s)),
            -s,
        )
    )
    defects = np.array(
        [
            [0.002, -0.001],
            [-0.003, 0.0015],
            [0.001, 0.002],
            [0.0, -0.001],
            [-0.002, -0.001],
            [0.0015, 0.001],
            [-0.001, 0.0025],
        ]
    )
    actual_ratio = ideal + defects

    ratio_error = float(np.max(np.linalg.norm(defects, axis=1)))
    first_moment = float(
        np.sum(weights * np.abs(v - midpoint) / length) / np.sum(weights)
    )
    cert = SignedPulseCovarianceErrorBudget(
        sigma=sigma,
        c_star=c_star,
        u_star=u_star,
        ratio_error_bound=math.nextafter(ratio_error, math.inf),
        normalized_first_moment_bound=math.nextafter(first_moment, math.inf),
    )

    normalized_column = _weighted_average(actual_ratio, weights)
    reference = np.asarray(cert.reference_direction)
    actual_error = normalized_column - reference

    assert float(np.linalg.norm(actual_error)) <= cert.normalized_error_bound
    assert cert.error_vector_for(normalized_column) == pytest.approx(actual_error)

    # Cross-check the analytic split independently from the production formula.
    average_defect = _weighted_average(defects, weights)
    center_drift = _weighted_average(ideal, weights) - reference
    assert actual_error == pytest.approx(average_defect + center_drift)
    assert float(np.linalg.norm(average_defect)) <= cert.ratio_error_bound
    assert float(np.linalg.norm(center_drift)) <= cert.concentration_error_upper


def test_global_ideal_direction_lipschitz_bound_covers_direct_derivative_norms():
    cert = SignedPulseCovarianceErrorBudget(
        sigma=-1,
        c_star=-1.25,
        u_star=2.0,
        ratio_error_bound=0.0,
        normalized_first_moment_bound=0.0,
    )

    for s in np.linspace(-20.0, 20.0, 101):
        derivative_norm = math.sqrt(
            1.0 + cert.c_star**2 * s**2 / (1.0 + s**2)
        )
        assert derivative_norm <= cert.ideal_direction_lipschitz_upper


def test_paper_component_rates_collapse_to_sqrt_scale_envelope():
    s_star = 625.0
    c_star = -0.8
    u_star = 1.6
    ratio_constant = 3.0
    first_moment_constant = 2.5

    ratio_bound = ratio_constant / s_star
    moment_bound = first_moment_constant / math.sqrt(s_star)
    cert = SignedPulseCovarianceErrorBudget(
        sigma=1,
        c_star=c_star,
        u_star=u_star,
        ratio_error_bound=ratio_bound,
        normalized_first_moment_bound=moment_bound,
    )
    envelope = sqrt_scale_error_envelope(
        s_star=s_star,
        c_star=c_star,
        u_star=u_star,
        ratio_constant=ratio_constant,
        first_moment_constant=first_moment_constant,
    )

    assert cert.normalized_error_bound <= envelope
    expected_constant = ratio_constant + u_star * math.hypot(1.0, c_star) * first_moment_constant
    assert envelope >= expected_constant / math.sqrt(s_star)


def test_uniform_two_sign_budget_is_ready_for_covariance_perturbation_delta():
    minus = SignedPulseCovarianceErrorBudget(
        sigma=-1,
        c_star=-0.6,
        u_star=1.3,
        ratio_error_bound=0.004,
        normalized_first_moment_bound=0.02,
    )
    plus = SignedPulseCovarianceErrorBudget(
        sigma=1,
        c_star=-0.6,
        u_star=1.3,
        ratio_error_bound=0.006,
        normalized_first_moment_bound=0.018,
    )

    assert uniform_two_sign_error_bound(minus, plus) == max(
        minus.normalized_error_bound,
        plus.normalized_error_bound,
    )

    with pytest.raises(ValueError, match="expected minus.sigma"):
        uniform_two_sign_error_bound(plus, minus)

    mismatched = SignedPulseCovarianceErrorBudget(
        sigma=1,
        c_star=-0.6,
        u_star=1.4,
        ratio_error_bound=0.006,
        normalized_first_moment_bound=0.018,
    )
    with pytest.raises(ValueError, match="identical c_star and u_star"):
        uniform_two_sign_error_bound(minus, mismatched)


def test_fail_closed_inputs_and_column_budget():
    with pytest.raises(ValueError, match="sigma"):
        SignedPulseCovarianceErrorBudget(
            sigma=0,
            c_star=-0.5,
            u_star=1.0,
            ratio_error_bound=0.0,
            normalized_first_moment_bound=0.0,
        )
    with pytest.raises(ValueError, match="c_star<0"):
        SignedPulseCovarianceErrorBudget(
            sigma=1,
            c_star=0.1,
            u_star=1.0,
            ratio_error_bound=0.0,
            normalized_first_moment_bound=0.0,
        )
    with pytest.raises(ValueError, match="nonnegative"):
        SignedPulseCovarianceErrorBudget(
            sigma=1,
            c_star=-0.5,
            u_star=1.0,
            ratio_error_bound=-1e-4,
            normalized_first_moment_bound=0.0,
        )

    cert = SignedPulseCovarianceErrorBudget(
        sigma=1,
        c_star=-0.5,
        u_star=1.0,
        ratio_error_bound=0.01,
        normalized_first_moment_bound=0.01,
    )
    with pytest.raises(ValueError, match="exceeds analytic"):
        cert.error_vector_for((10.0, 10.0))
    with pytest.raises(ValueError, match="s_star>=1"):
        sqrt_scale_error_envelope(
            s_star=0.5,
            c_star=-0.5,
            u_star=1.0,
            ratio_constant=1.0,
            first_moment_constant=1.0,
        )
