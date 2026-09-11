import math

import numpy as np
import pytest

from openai_ns_reconstruction.section10_support_energy import (
    FIXED_TIME_ENERGY_COEFFICIENT_OVER_PI_EXACT,
    SUPPORT_HALF_HEIGHT_EXACT,
    SUPPORT_HEIGHT_EXACT,
    SUPPORT_RADIUS_EXACT,
    SUPPORT_RADIUS_SQUARED_EXACT,
    SUPPORT_VOLUME_OVER_PI_EXACT,
    Section10SupportEnergyCertificate,
    fixed_time_energy_bound_from_speed_sup,
    geometry_matches_spatial_localization_constants,
    section10_support_volume,
)
from openai_ns_reconstruction.spatial_localization import section10_axisymmetric_cutoff
from openai_ns_reconstruction.verify import kinetic_energy_axisymmetric


def test_exact_support_geometry_and_energy_coefficient():
    assert SUPPORT_RADIUS_SQUARED_EXACT.numerator == 1
    assert SUPPORT_RADIUS_SQUARED_EXACT.denominator == 16
    assert SUPPORT_RADIUS_EXACT * SUPPORT_RADIUS_EXACT == SUPPORT_RADIUS_SQUARED_EXACT
    assert SUPPORT_HALF_HEIGHT_EXACT * 2 == SUPPORT_HEIGHT_EXACT
    assert SUPPORT_VOLUME_OVER_PI_EXACT == SUPPORT_RADIUS_SQUARED_EXACT * SUPPORT_HEIGHT_EXACT
    assert FIXED_TIME_ENERGY_COEFFICIENT_OVER_PI_EXACT == SUPPORT_VOLUME_OVER_PI_EXACT / 2
    assert section10_support_volume() == pytest.approx(math.pi / 32.0)
    assert geometry_matches_spatial_localization_constants()


def test_fixed_time_energy_bound_is_exact_geometry_propagation():
    cert = Section10SupportEnergyCertificate(3.0)
    assert cert.support_radius == 0.25
    assert cert.support_half_height == 0.25
    assert cert.support_volume == pytest.approx(math.pi / 32.0)
    assert cert.fixed_time_energy_upper_bound == pytest.approx(9.0 * math.pi / 64.0)
    assert fixed_time_energy_bound_from_speed_sup(0.0) == 0.0

    assert cert.enclosing_support_contains(0.0, 0.0, 0.0)
    assert cert.enclosing_support_contains(0.25, 0.0, 0.0)
    assert not cert.enclosing_support_contains(0.250001, 0.0, 0.0)
    assert not cert.enclosing_support_contains(0.0, 0.0, 0.250001)


def test_positive_quadrature_fixture_respects_analytic_energy_bound():
    # Independent test fixture only: it checks the support-volume implication,
    # not the unresolved paper velocity.  Since 0<=c<=1, |u|<=M pointwise.
    M = 2.75

    def bounded_fixture(r: float, z: float, t: float) -> np.ndarray:
        c = section10_axisymmetric_cutoff(r, z, t)
        return np.array([0.0, 0.0, M * c])

    energy = kinetic_energy_axisymmetric(
        bounded_fixture,
        0.4,
        radius=0.25,
        z_min=-0.25,
        z_max=0.25,
        n=24,
    )
    bound = fixed_time_energy_bound_from_speed_sup(M)
    assert 0.0 < energy <= bound * (1.0 + 1.0e-13)


def test_certificate_rejects_uncertified_nonfinite_or_negative_bounds():
    for bad in (-1.0, float("nan"), float("inf"), float("-inf")):
        with pytest.raises(ValueError):
            Section10SupportEnergyCertificate(bad)
        with pytest.raises(ValueError):
            fixed_time_energy_bound_from_speed_sup(bad)
