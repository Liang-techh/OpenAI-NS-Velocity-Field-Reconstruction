import math

import numpy as np
import pytest

from openai_ns_reconstruction.local_field import LocalField, LocalizedField
from openai_ns_reconstruction.spatial_localization import (
    PLATEAU_HALF_HEIGHT,
    PLATEAU_RADIUS_SQUARED,
    SUPPORT_HALF_HEIGHT,
    SUPPORT_RADIUS_SQUARED,
    plateau_contains,
    section10_axisymmetric_cutoff,
    section10_spatial_cutoff,
    support_cylinder_contains,
    support_is_strictly_inside_period_cube,
    symmetric_smooth_bump,
)


def test_symmetric_bump_matches_formal_support_and_plateau_axioms():
    for s in (-2.0, -1.0, -0.75, -0.5, -0.1, 0.0, 0.1, 0.5, 0.75, 1.0, 2.0):
        assert symmetric_smooth_bump(s) == pytest.approx(symmetric_smooth_bump(-s))
    assert symmetric_smooth_bump(-0.5) == 1.0
    assert symmetric_smooth_bump(0.5) == 1.0
    assert symmetric_smooth_bump(-1.0) == 0.0
    assert symmetric_smooth_bump(1.0) == 0.0


def test_section10_cutoff_has_official_plateau_and_closed_support_geometry():
    r_plateau = math.sqrt(PLATEAU_RADIUS_SQUARED) * 0.99
    assert plateau_contains(r_plateau, 0.0, PLATEAU_HALF_HEIGHT * 0.99)
    assert section10_spatial_cutoff(r_plateau, 0.0, PLATEAU_HALF_HEIGHT * 0.99) == 1.0

    # Boundary and exterior points are exactly zero for the executable bump.
    r_support = math.sqrt(SUPPORT_RADIUS_SQUARED)
    assert support_cylinder_contains(r_support, 0.0, 0.0)
    assert section10_spatial_cutoff(r_support, 0.0, 0.0) == 0.0
    assert support_cylinder_contains(0.0, 0.0, SUPPORT_HALF_HEIGHT)
    assert section10_spatial_cutoff(0.0, 0.0, SUPPORT_HALF_HEIGHT) == 0.0
    assert section10_spatial_cutoff(0.0, 0.0, -SUPPORT_HALF_HEIGHT) == 0.0
    assert not support_cylinder_contains(r_support + 1e-6, 0.0, 0.0)
    assert not support_cylinder_contains(0.0, 0.0, SUPPORT_HALF_HEIGHT + 1e-6)


def test_section10_cutoff_is_axisymmetric_and_inside_unit_period_cube():
    r, z, t = 0.13, 0.08, 0.6
    reference = section10_axisymmetric_cutoff(r, z, t)
    for theta in np.linspace(0.0, 2.0 * math.pi, 17):
        x, y = r * math.cos(theta), r * math.sin(theta)
        assert section10_spatial_cutoff(x, y, z, t) == pytest.approx(reference)
        assert support_is_strictly_inside_period_cube(x, y, z)

    assert not support_is_strictly_inside_period_cube(0.3, 0.0, 0.0)


def test_negative_z_exterior_is_zero_and_short_circuits_missing_local_field():
    # Regression guard: the old one-sided smooth_cutoff(s) is not suitable by
    # itself for the official |z|-symmetric compact support.
    assert section10_spatial_cutoff(0.0, 0.0, -0.30, 0.5) == 0.0

    def unavailable(*_args):
        raise RuntimeError("local field must not be evaluated outside support")

    local = LocalField(unavailable, unavailable)
    localized = LocalizedField.from_axisymmetric(local, section10_axisymmetric_cutoff)
    assert np.array_equal(localized.velocity(0.0, 0.0, -0.30, 0.5), np.zeros(3))
    assert np.array_equal(localized.velocity(0.30, 0.0, 0.0, 0.5), np.zeros(3))


def test_bad_geometry_inputs_fail_loudly():
    with pytest.raises(ValueError):
        section10_spatial_cutoff(float("nan"), 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        section10_axisymmetric_cutoff(-0.1, 0.0, 0.0)
