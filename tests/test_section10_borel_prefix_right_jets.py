import math
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_jets import (
    borel_extension_from_full_spacetime_jets,
)
from openai_ns_reconstruction.section10_borel_prefix_right_jets import (
    certify_section10_borel_prefix_right_jets,
)


NORMAL_JETS = (
    np.array([1.2, -0.4, 0.7]),
    np.array([-0.3, 0.8, 0.2]),
    np.array([0.6, 0.1, -0.5]),
    np.array([-0.2, 0.35, 0.45]),
)


def _finite_support_full_jets(degree: int, x: float, y: float, z: float) -> np.ndarray:
    shape = (3,) + (4,) * degree
    out = np.zeros(shape)
    value = NORMAL_JETS[degree] if degree < len(NORMAL_JETS) else np.zeros(3)
    # Keep a harmless spatial dependence only in degree zero so the test also
    # exercises evaluation at the requested point without changing time jets.
    if degree == 0:
        value = value + np.array([0.1 * x, -0.2 * y, 0.05 * z])
    out[(slice(None),) + (0,) * degree] = value
    return out


def test_prefix_certificate_matches_independent_future_branch_interpolation() -> None:
    """Recover endpoint jets by sampling the actual landed Borel evaluator.

    Production uses exact Taylor-prefix algebra.  This regression instead
    evaluates the future branch at four positive offsets inside the common
    plateau, solves a Vandermonde interpolation problem, and differentiates the
    independently recovered cubic polynomial at the endpoint.
    """

    point = (0.17, -0.09, 0.23)
    local_scale = lambda degree: 1 << degree
    certificate = certify_section10_borel_prefix_right_jets(
        _finite_support_full_jets,
        local_scale,
        *point,
        3,
    )
    assert certificate.formal_prefix_ready

    right = borel_extension_from_full_spacetime_jets(
        _finite_support_full_jets,
        local_scale,
        endpoint=1.0,
    )

    # The degree-three common plateau radius is 1/(2*8)=1/16.  All four
    # interpolation points stay well inside it.  Higher endpoint coefficients
    # are exactly zero in this analytic fixture, so the evaluated future branch
    # is a cubic on this neighborhood.
    radius = float(certificate.rows[3].common_plateau_radius)
    step = radius / 8.0
    offsets = np.arange(4, dtype=float) * step
    values = np.vstack(
        [right(*point, 1.0 + float(offset)) for offset in offsets]
    )
    vandermonde = np.vander(offsets, N=4, increasing=True)
    coefficients = np.linalg.solve(vandermonde, values)

    for degree, row in enumerate(certificate.rows):
        recovered = math.factorial(degree) * coefficients[degree]
        assert recovered == pytest.approx(
            row.prefix_endpoint_time_jet,
            rel=2e-9,
            abs=2e-9,
        )
        assert row.prefix_endpoint_time_jet == pytest.approx(
            row.prescribed_normal_jet,
            rel=0.0,
            abs=0.0,
        )
        assert row.finite_prefix_right_jet_verified


def test_common_plateau_radius_uses_exact_doubling_envelope_arithmetic() -> None:
    certificate = certify_section10_borel_prefix_right_jets(
        _finite_support_full_jets,
        lambda degree: 3 * (1 << degree),
        0.0,
        0.0,
        0.0,
        3,
    )

    assert [row.prefix_scales for row in certificate.rows] == [
        (3,),
        (3, 6),
        (3, 6, 12),
        (3, 6, 12, 24),
    ]
    assert [row.common_plateau_radius for row in certificate.rows] == [
        Fraction(1, 6),
        Fraction(1, 12),
        Fraction(1, 24),
        Fraction(1, 48),
    ]
    for row in certificate.rows:
        assert row.unit_plateau_verified
        for scale in row.prefix_scales:
            assert Fraction(scale, 1) * row.common_plateau_radius <= Fraction(1, 2)


def test_wrong_full_jet_shape_and_non_section10_endpoint_fail_closed() -> None:
    def wrong_shape(degree: int, x: float, y: float, z: float) -> np.ndarray:
        if degree == 2:
            return np.zeros((3, 4))
        return _finite_support_full_jets(degree, x, y, z)

    with pytest.raises(ValueError, match="order-2 full spacetime jet"):
        certify_section10_borel_prefix_right_jets(
            wrong_shape,
            lambda degree: 1 << degree,
            0.0,
            0.0,
            0.0,
            2,
        )

    with pytest.raises(ValueError, match="requires endpoint=1"):
        certify_section10_borel_prefix_right_jets(
            _finite_support_full_jets,
            lambda degree: 1 << degree,
            0.0,
            0.0,
            0.0,
            2,
            endpoint=0.9,
        )


def test_certificate_keeps_infinite_smoothness_and_force_boundaries_closed() -> None:
    certificate = certify_section10_borel_prefix_right_jets(
        _finite_support_full_jets,
        lambda degree: 1 << degree,
        0.0,
        0.0,
        0.0,
        2,
    )

    assert certificate.dense_full_jet_shape_verified
    assert certificate.finite_prefix_right_jet_algebra_verified
    assert not certificate.actual_section9_residual_limits_verified
    assert not certificate.analytic_template_bounds_verified
    assert not certificate.infinite_borel_right_jets_verified
    assert not certificate.all_order_borel_smoothness_verified
    assert not certificate.smooth_compact_forcing_verified
    assert not certificate.paper_exact_velocity_available


def test_inputs_reject_lossy_or_invalid_degree_and_scale_contracts() -> None:
    for bad_degree in (-1, True, 1.5):
        with pytest.raises(ValueError, match="nonnegative integer"):
            certify_section10_borel_prefix_right_jets(
                _finite_support_full_jets,
                lambda degree: 1 << degree,
                0.0,
                0.0,
                0.0,
                bad_degree,
            )

    with pytest.raises(ValueError, match="local_scale must be callable"):
        certify_section10_borel_prefix_right_jets(
            _finite_support_full_jets,
            None,
            0.0,
            0.0,
            0.0,
            1,
        )
