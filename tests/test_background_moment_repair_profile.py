import math

import numpy as np
import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


def _moment_jets(eta):
    eta = float(eta)
    values = np.array(
        [
            0.20 + 0.10 * eta,
            -0.15 + 0.02 * eta**2,
            0.11 - 0.03 * eta,
            0.25 + 0.04 * eta + 0.01 * eta**2,
            -0.12 + 0.05 * eta,
        ]
    )
    first = np.array(
        [
            0.10,
            0.04 * eta,
            -0.03,
            0.04 + 0.02 * eta,
            0.05,
        ]
    )
    return np.vstack((values, first))


def _patch_jets(eta):
    return np.array([1.5 + 0.2 * float(eta), 0.2])


def _adapter():
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )


def _integrate_on_bumps(bumps, integrand):
    return sum(
        quad(integrand, bump.left, bump.right, epsabs=2e-12, epsrel=2e-12)[0]
        for bump in bumps
    )


def test_functional_repair_cancels_supplied_moments_by_direct_profile_quadrature():
    adapter = _adapter()
    repair = adapter.repair
    base = adapter.base_profile
    eta = 0.31
    patch = _patch_jets(eta)[0]
    lam = repair.lambda_exponent

    def u_corr(radius):
        X = 0.5 * radius * radius
        return adapter.U(X, eta) - base.U(X, eta)

    def e_corr(radius):
        X = 0.5 * radius * radius
        return adapter.E(X, eta) - base.E(X, eta)

    correction = np.array(
        [
            _integrate_on_bumps(repair.u_bumps, lambda r: r * u_corr(r)),
            _integrate_on_bumps(repair.e_bumps, lambda r: r**2 * e_corr(r)),
            2.0
            * patch
            * _integrate_on_bumps(
                repair.e_bumps,
                lambda r: r ** (-2.0 - 2.0 * lam) * e_corr(r),
            ),
            patch
            * _integrate_on_bumps(
                repair.u_bumps,
                lambda r: r ** (1.0 - 2.0 * lam) * u_corr(r),
            ),
            -patch
            * _integrate_on_bumps(
                repair.e_bumps,
                lambda r: r ** (-2.0 * lam) * e_corr(r),
            ),
        ]
    )
    assert np.allclose(_moment_jets(eta)[0] + correction, 0.0, rtol=0.0, atol=3e-10)


def test_repaired_dU_deta_matches_independent_parameter_difference():
    adapter = _adapter()
    eta = -0.27
    X = 0.5 * 2.17**2
    step = 2.0e-5
    finite_difference = (
        adapter.U(X, eta + step) - adapter.U(X, eta - step)
    ) / (2.0 * step)
    assert math.isclose(
        adapter.dU_deta(X, eta),
        finite_difference,
        rel_tol=5e-6,
        abs_tol=2e-8,
    )


def test_compact_functional_repair_preserves_base_profile_off_patch_and_gate():
    adapter = _adapter()
    base = adapter.base_profile
    eta = 0.2
    for radius in (0.0, 1.0, 7.0):
        X = 0.5 * radius * radius
        assert adapter.U(X, eta) == base.U(X, eta)
        assert adapter.E(X, eta) == base.E(X, eta)
        assert adapter.dU_deta(X, eta) == base.dU_deta(X, eta)

    profile = adapter.as_leading_profile()
    assert profile.paper_exact is False
    assert profile.Pi is None
    assert "formal-structure" in profile.provenance
    assert profile.smooth_swirl_factor(0.0, eta) == base.smooth_swirl_factor(0.0, eta)


def test_functional_repair_fails_closed_on_missing_derivative_or_zero_patch_factor():
    base = toy_gaussian_profile()
    repair = _repair()
    missing_derivative = Lemma52RepairedProfileAdapter(
        repair,
        base,
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    with pytest.raises(ValueError):
        missing_derivative.dU_deta(0.5 * 2.1**2, 0.1)

    zero_patch = Lemma52RepairedProfileAdapter(
        repair,
        base,
        _moment_jets,
        lambda eta: np.array([0.0, 0.2]),
    )
    with pytest.raises(ValueError):
        zero_patch.U(0.5 * 2.1**2, 0.1)
