import math

import numpy as np
import pytest

from openai_ns_reconstruction.section10_origin_preservation import (
    BLOWUP_TIME,
    ORIGIN,
    Section10LateOriginCertificate,
    certify_section10_late_origin,
    late_origin_cut_then_time_value,
)
from openai_ns_reconstruction.spatial_localization import (
    plateau_contains,
    section10_spatial_cutoff,
    section10_spatial_cutoff_gradient,
)
from openai_ns_reconstruction.time_localization import LATE_START, time_switch


@pytest.mark.parametrize(
    "t",
    [LATE_START, 0.9, np.nextafter(BLOWUP_TIME, 0.0)],
)
def test_late_origin_certificate_matches_independent_localization_factors(t):
    cert = certify_section10_late_origin(t)

    # Cross-module oracle: inspect the already-landed spatial/time adapters
    # directly rather than trusting only fields cached by the certificate.
    assert plateau_contains(*ORIGIN)
    assert section10_spatial_cutoff(*ORIGIN, t) == 1.0
    assert np.array_equal(section10_spatial_cutoff_gradient(*ORIGIN, t), np.zeros(3))
    assert time_switch(t) == 1.0

    assert cert.origin_in_plateau
    assert cert.spatial_cutoff_value == 1.0
    assert np.array_equal(cert.spatial_cutoff_gradient, np.zeros(3))
    assert cert.time_switch_value == 1.0
    assert cert.identity_factors_exact


def test_cut_product_rule_and_time_activation_preserve_arbitrary_finite_origin_data():
    # Algebra fixture only: these vectors are deliberately not represented as
    # a paper field.  The theorem-shaped implication is that *whatever* the
    # upstream certified A and curl(A) values are, late localization does not
    # alter curl(A) at the origin.
    curl_value = np.array([2.0, -3.0, 5.0])
    potential_value = np.array([7.0, 11.0, -13.0])
    cert = Section10LateOriginCertificate(0.875)

    cut_value = cert.cut_velocity_from_potential_data(curl_value, potential_value)
    activated = cert.activate_periodic_velocity(cut_value)
    composed = late_origin_cut_then_time_value(curl_value, potential_value, 0.875)

    assert np.array_equal(cut_value, curl_value)
    assert np.array_equal(activated, curl_value)
    assert np.array_equal(composed, curl_value)


def test_certificate_fails_closed_outside_left_limit_late_branch():
    invalid = [
        np.nextafter(LATE_START, 0.0),
        BLOWUP_TIME,
        np.nextafter(BLOWUP_TIME, math.inf),
        float("nan"),
        float("inf"),
        float("-inf"),
    ]
    for t in invalid:
        with pytest.raises(ValueError):
            Section10LateOriginCertificate(t)


def test_point_data_must_be_finite_length_three_vectors():
    cert = Section10LateOriginCertificate(0.8)

    with pytest.raises(ValueError):
        cert.cut_velocity_from_potential_data([1.0, 2.0], [0.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        cert.cut_velocity_from_potential_data(
            [1.0, 2.0, float("nan")], [0.0, 0.0, 0.0]
        )
    with pytest.raises(ValueError):
        cert.activate_periodic_velocity([0.0, float("inf"), 0.0])
