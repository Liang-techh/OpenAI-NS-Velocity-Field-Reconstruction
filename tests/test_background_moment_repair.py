import math

import numpy as np
import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.background_moment_repair import (
    CompactMomentBump,
    Lemma52MomentRepair,
)


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


def test_compact_bumps_are_unit_mass_nonnegative_and_flatly_supported():
    repair = _repair()
    for bump in repair.u_bumps + repair.e_bumps:
        mass = quad(lambda r: bump(r), bump.left, bump.right, epsabs=1e-13)[0]
        assert math.isclose(mass, 1.0, rel_tol=2e-12, abs_tol=2e-12)
        assert bump(bump.left) == 0.0
        assert bump(bump.right) == 0.0
        assert bump(0.5 * (bump.left + bump.right)) > 0.0
        assert bump(bump.left - 0.1) == 0.0
        assert bump(bump.right + 0.1) == 0.0


def test_eq_5_16_cancels_all_five_moments_with_independent_quadrature():
    repair = _repair()
    base = np.array([0.31, -0.27, 0.19, 0.41, -0.23])
    patch = 1.7
    out = repair.solve(base, patch)

    def u_corr(r):
        return sum(a * b(r) for a, b in zip(out.alpha, repair.u_bumps))

    def e_corr(r):
        return sum(bcoef * bump(r) for bcoef, bump in zip(out.beta, repair.e_bumps))

    u_left = min(b.left for b in repair.u_bumps)
    u_right = max(b.right for b in repair.u_bumps)
    e_left = min(b.left for b in repair.e_bumps)
    e_right = max(b.right for b in repair.e_bumps)
    lam = repair.lambda_exponent

    # Independent adaptive quadrature: production assembles B_U/B_E with the
    # repository's cached Gauss-Legendre rule.
    correction = np.array([
        quad(lambda r: r * u_corr(r), u_left, u_right, epsabs=1e-12)[0],
        quad(lambda r: r**2 * e_corr(r), e_left, e_right, epsabs=1e-12)[0],
        2.0 * patch * quad(
            lambda r: r**(-2.0 - 2.0 * lam) * e_corr(r),
            e_left,
            e_right,
            epsabs=1e-12,
        )[0],
        patch * quad(
            lambda r: r**(1.0 - 2.0 * lam) * u_corr(r),
            u_left,
            u_right,
            epsabs=1e-12,
        )[0],
        -patch * quad(
            lambda r: r**(-2.0 * lam) * e_corr(r),
            e_left,
            e_right,
            epsabs=1e-12,
        )[0],
    ])
    assert np.allclose(base + correction, 0.0, rtol=0.0, atol=2e-10)
    assert np.max(np.abs(out.corrected_moments)) < 1e-11


def test_eq_5_14_correction_is_compact_and_leaves_inner_axis_values_unchanged():
    repair = _repair()
    out = repair.solve([0.2, -0.1, 0.05, 0.3, -0.2], 2.0)
    base_u, base_e = 3.25, -1.5

    assert repair.corrected_u_value(0.0, base_u, out) == base_u
    assert repair.corrected_e_value(0.0, base_e, out) == base_e
    assert repair.corrected_u_value(0.5, base_u, out) == base_u
    assert repair.corrected_e_value(0.5, base_e, out) == base_e

    outer_X = 0.5 * 7.0**2
    assert repair.corrected_u_value(outer_X, base_u, out) == base_u
    assert repair.corrected_e_value(outer_X, base_e, out) == base_e

    u_X = 0.5 * 2.15**2
    e_X = 0.5 * 4.15**2
    assert repair.corrected_u_value(u_X, base_u, out) != base_u
    assert repair.corrected_e_value(e_X, base_e, out) != base_e


def test_moment_matrices_use_the_paper_exponents_and_are_invertible():
    repair = _repair()
    lam = repair.lambda_exponent
    expected_u = np.array([
        [quad(lambda r, b=b: r * b(r), b.left, b.right)[0] for b in repair.u_bumps],
        [
            quad(lambda r, b=b: r**(1 - 2 * lam) * b(r), b.left, b.right)[0]
            for b in repair.u_bumps
        ],
    ])
    expected_e = np.array([
        [quad(lambda r, b=b: r**2 * b(r), b.left, b.right)[0] for b in repair.e_bumps],
        [
            quad(lambda r, b=b: r**(-2 - 2 * lam) * b(r), b.left, b.right)[0]
            for b in repair.e_bumps
        ],
        [
            quad(lambda r, b=b: r**(-2 * lam) * b(r), b.left, b.right)[0]
            for b in repair.e_bumps
        ],
    ])
    assert np.allclose(repair.u_matrix, expected_u, rtol=2e-11, atol=2e-12)
    assert np.allclose(repair.e_matrix, expected_e, rtol=2e-11, atol=2e-12)
    assert abs(np.linalg.det(repair.u_matrix)) > 1e-3
    assert abs(np.linalg.det(repair.e_matrix)) > 1e-4


@pytest.mark.parametrize(
    "factory",
    [
        lambda: Lemma52MomentRepair.from_intervals(
            0.0, ((1, 2), (3, 4)), ((5, 6), (7, 8), (9, 10))
        ),
        lambda: Lemma52MomentRepair.from_intervals(
            0.2, ((1, 3), (2, 4)), ((5, 6), (7, 8), (9, 10))
        ),
        lambda: Lemma52MomentRepair.from_intervals(
            0.2, ((1, 2),), ((5, 6), (7, 8), (9, 10))
        ),
    ],
)
def test_bad_repair_geometry_fails_closed(factory):
    with pytest.raises((ValueError, TypeError)):
        factory()


@pytest.mark.parametrize(
    "moments,patch",
    [
        ([1, 2, 3, 4], 1.0),
        ([1, 2, 3, 4, math.nan], 1.0),
        ([1, 2, 3, 4, 5], 0.0),
        ([1, 2, 3, 4, 5], math.inf),
    ],
)
def test_bad_moment_data_fails_closed(moments, patch):
    with pytest.raises(ValueError):
        _repair().solve(moments, patch)


def test_bad_profile_coordinate_fails_closed():
    repair = _repair()
    out = repair.solve([1, 2, 3, 4, 5], 2.0)
    with pytest.raises(ValueError):
        repair.corrected_u_value(-1.0, 0.0, out)
    with pytest.raises(ValueError):
        repair.corrected_e_value(math.nan, 0.0, out)
