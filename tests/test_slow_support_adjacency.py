import math

import pytest

from openai_ns_reconstruction.charts import DyadicChart
from openai_ns_reconstruction.slow_labels import SlowLabel
from openai_ns_reconstruction.slow_support_adjacency import (
    PhysicalSlowBox,
    axis_exponents,
    certify_slow_support_adjacency,
    physical_grid_center,
    physical_mesh_widths,
)


def _nearest_grid(chart, R, Z, T):
    mesh = chart.S_star ** -3
    return tuple(int(math.floor(value / mesh + 0.5)) for value in (R, Z, T))


def _interior_slow_point(chart, *, s=0.8, R=1.0, Z=0.1):
    # Normalized form of q-z^2 q^(2h)=tau with q=Q*s:
    # s - Z^2 s^(2h) = T.
    T = s - Z * Z * s ** (2.0 * chart.h)
    assert T > 0.0
    return R, Z, T


def test_physical_mesh_widths_match_pinned_slotcoloring_spacing_formula():
    h = 0.005
    label = SlowLabel(ell=100, a=(10, -3, 7), sigma=1)
    widths = physical_mesh_widths(label, h=h)

    for width, exponent in zip(widths, axis_exponents(h)):
        independent = 2.0 ** (-label.ell * exponent) / label.ell ** 6
        assert width == pytest.approx(independent, rel=8e-15, abs=0.0)

    center = physical_grid_center(label, h=h)
    assert center == pytest.approx(tuple(w * a for w, a in zip(widths, label.a)), rel=0.0, abs=0.0)


def test_landed_one_mesh_product_support_is_inside_official_two_mesh_box():
    h = 0.005
    chart = DyadicChart(100, h)
    mesh = chart.S_star ** -3
    label = SlowLabel(
        ell=chart.ell,
        a=_nearest_grid(chart, 1.0, 0.2, 0.8),
        sigma=1,
    )
    box = PhysicalSlowBox(label, h)
    center = tuple(mesh * a for a in label.a)

    normalized_inside = (
        center[0] + 0.9 * mesh,
        center[1] - 0.4 * mesh,
        center[2] + 0.75 * mesh,
    )
    physical_inside = chart.to_physical_tau(*normalized_inside)
    assert box.contains_landed_cutoff_support_point(*physical_inside)
    assert box.contains(*physical_inside)

    normalized_outside = (center[0] + 2.1 * mesh, center[1], center[2])
    physical_outside = chart.to_physical_tau(*normalized_outside)
    assert not box.contains(*physical_outside)


def test_cross_band_support_overlap_uses_same_physical_point_and_certifies_slots():
    h = 0.005
    left_chart = DyadicChart(100, h)
    R, Z, T = _interior_slow_point(left_chart, s=0.8)
    physical = left_chart.to_physical_tau(R, Z, T)

    left = SlowLabel(left_chart.ell, _nearest_grid(left_chart, R, Z, T), 1)
    right_chart = DyadicChart(101, h)
    R1, Z1, T1 = right_chart.from_physical_tau(*physical)
    right = SlowLabel(right_chart.ell, _nearest_grid(right_chart, R1, Z1, T1), -1)

    cert = certify_slow_support_adjacency(
        left, right, r=physical[0], z=physical[1], tau=physical[2], h=h
    )
    assert cert.certified
    assert cert.level_gap == 1
    assert cert.slot.certified
    assert cert.q / left_chart.Q == pytest.approx(0.8, rel=3e-13)
    assert 0.5 < cert.left_q_ratio < 2.0
    assert 0.5 < cert.right_q_ratio < 2.0
    assert all(abs(x) <= 0.5 + 2e-4 for x in cert.left_offsets)
    assert all(abs(x) <= 0.5 + 2e-4 for x in cert.right_offsets)


def test_same_level_two_mesh_enlargement_recovers_grid_gap_constraint():
    h = 0.005
    chart = DyadicChart(100, h)
    mesh = chart.S_star ** -3
    _, Z, T = _interior_slow_point(chart, s=0.8)
    radial_mid = int(math.floor(1.0 / mesh + 0.5))
    z_idx = int(math.floor(Z / mesh + 0.5))
    t_idx = int(math.floor(T / mesh + 0.5))

    # Use an interior common point rather than the exact closed-box boundary:
    # binary64 round trips are not interval certificates.  The exact gap-four
    # discrete boundary is already tested independently in test_slot_geometry.
    left = SlowLabel(chart.ell, (radial_mid - 1, z_idx, t_idx), 1)
    right = SlowLabel(chart.ell, (radial_mid + 2, z_idx, t_idx), 1)
    normalized = ((radial_mid + 0.25) * mesh, z_idx * mesh, t_idx * mesh)
    physical = chart.to_physical_tau(*normalized)

    cert = certify_slow_support_adjacency(
        left,
        right,
        r=physical[0],
        z=physical[1],
        tau=physical[2],
        h=h,
        product_radius_meshes=2.0,
    )
    assert cert.certified
    assert abs(left.a[0] - right.a[0]) == 3
    assert cert.slot.lower_color != cert.slot.upper_color

    with pytest.raises(ValueError, match="product-support enclosure"):
        certify_slow_support_adjacency(
            left,
            right,
            r=physical[0],
            z=physical[1],
            tau=physical[2],
            h=h,
            product_radius_meshes=1.0,
        )


def test_adjacency_bridge_fails_closed_outside_landed_support_hypotheses():
    h = 0.005
    chart = DyadicChart(100, h)
    R, Z, T = _interior_slow_point(chart, s=0.8)
    physical = chart.to_physical_tau(R, Z, T)
    left = SlowLabel(chart.ell, _nearest_grid(chart, R, Z, T), 1)

    far_chart = DyadicChart(103, h)
    R3, Z3, T3 = far_chart.from_physical_tau(*physical)
    far = SlowLabel(far_chart.ell, _nearest_grid(far_chart, R3, Z3, T3), -1)
    with pytest.raises(ValueError, match="dyadic cutoff support"):
        certify_slow_support_adjacency(
            left, far, r=physical[0], z=physical[1], tau=physical[2], h=h
        )

    with pytest.raises(ValueError, match=r"\(0,2\]"):
        certify_slow_support_adjacency(
            left,
            left.opposite(),
            r=physical[0],
            z=physical[1],
            tau=physical[2],
            h=h,
            product_radius_meshes=2.01,
        )
