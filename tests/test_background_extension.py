import math

import numpy as np
import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.background_extension import reconstruct_eq_5_15
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair


def test_eq_5_15_matches_independent_direct_integrals():
    h = 0.005
    n = 1
    C = 2.75
    X = 0.8
    eta = 0.3

    def U(x, e):
        return 1.0 + e + (2.0 - e) * x + 0.5 * x * x

    def dU_deta(x, e):
        return 1.0 - x

    def dU_dX(x, e):
        return 2.0 - e + x

    phi0 = lambda x, e: 1.2 + 0.3 * x + 0.1 * e
    phi1 = lambda x, e: -0.4 + 0.2 * x - 0.05 * e
    omega_over_x = lambda x, e: 0.6 - 0.1 * x + 0.02 * e

    out = reconstruct_eq_5_15(
        n,
        X,
        eta,
        U,
        dU_deta,
        (phi0, phi1),
        omega_over_x,
        h=h,
        C=C,
        quadrature_points=64,
    )

    # Independent SciPy quadrature of the three formulas printed in Eq. (5.15).
    expected_F = quad(lambda x: U(x, eta), 0.0, X, epsabs=1e-13)[0]
    expected_dF = quad(lambda x: dU_deta(x, eta), 0.0, X, epsabs=1e-13)[0]

    A = 0.5 + h
    lam = 2.0 * n * h
    b = -A + lam
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta

    def direct_Z(x):
        return (
            2.0 * b * eta * U(x, eta)
            + d * dU_deta(x, eta)
            - 2.0 * eta * x * dU_dX(x, eta)
        ) / L

    expected_V = -quad(direct_Z, 0.0, X, epsabs=1e-13)[0]
    expected_Pi = quad(
        lambda x: (
            2.0 * phi0(x, eta) * phi1(x, eta) / (C * C)
            - 0.5 * omega_over_x(x, eta)
        ),
        0.0,
        X,
        epsabs=1e-13,
    )[0]

    assert math.isclose(out.F, expected_F, rel_tol=2e-13, abs_tol=2e-13)
    assert math.isclose(out.dF_deta, expected_dF, rel_tol=2e-13, abs_tol=2e-13)
    assert math.isclose(out.V, expected_V, rel_tol=2e-12, abs_tol=2e-12)
    assert math.isclose(out.Pi, expected_Pi, rel_tol=2e-12, abs_tol=2e-12)


def test_eq_5_15_consumes_actual_compact_repair_outputs():
    repair = Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )
    coefficients = repair.solve([0.31, -0.27, 0.19, 0.41, -0.23], 1.7)
    C = 3.0
    eta = 0.25
    X = 0.5 * 3.2**2

    def base_U(x, e):
        return 0.4 + 0.2 * e + 0.03 * x

    def U(x, e):
        return repair.corrected_u_value(x, base_U(x, e), coefficients)

    # The selected repair data are constant in eta in this test, so the compact
    # correction has zero eta derivative.  This is an interoperability check,
    # not a claim about the paper's still-missing coefficient functions.
    def dU_deta(x, e):
        return 0.2

    def base_phi1(x, e):
        return 0.35 + 0.01 * x

    def phi1(x, e):
        if x == 0.0:
            return base_phi1(0.0, e)
        radius = math.sqrt(2.0 * x)
        base_E = radius * base_phi1(x, e) / C
        repaired_E = repair.corrected_e_value(x, base_E, coefficients)
        return C * repaired_E / radius

    phi0 = lambda x, e: 1.0 + 0.02 * x
    omega_over_x = lambda x, e: 0.15 + 0.01 * x

    # The compact C-infinity repair bumps are deliberately flat at their joins.
    # Use a resolved fixed rule here and compare it to adaptive SciPy quadrature;
    # the earlier 128-node rule under-resolved the first narrow transition by
    # about 2.4e-5, which the independent oracle correctly exposed.
    out = reconstruct_eq_5_15(
        1,
        X,
        eta,
        U,
        dU_deta,
        (phi0, phi1),
        omega_over_x,
        h=0.005,
        C=C,
        quadrature_points=512,
    )

    expected_F = quad(lambda x: U(x, eta), 0.0, X, epsabs=2e-11, limit=200)[0]
    expected_Pi = quad(
        lambda x: 2.0 * phi0(x, eta) * phi1(x, eta) / (C * C)
        - 0.5 * omega_over_x(x, eta),
        0.0,
        X,
        epsabs=2e-10,
        limit=200,
    )[0]
    assert math.isclose(out.F, expected_F, rel_tol=2e-10, abs_tol=2e-10)
    assert math.isclose(out.Pi, expected_Pi, rel_tol=2e-10, abs_tol=2e-10)
    assert np.all(np.isfinite([out.F, out.V, out.Pi, out.dF_deta]))


def test_eq_5_15_axis_data_are_zero_by_the_forward_integrals():
    out = reconstruct_eq_5_15(
        1,
        0.0,
        -0.4,
        lambda x, e: 1.0 + e,
        lambda x, e: 1.0,
        (lambda x, e: 2.0, lambda x, e: 0.5),
        lambda x, e: 0.25,
        h=0.005,
        C=2.0,
    )
    assert out.F == 0.0
    assert out.V == 0.0
    assert out.Pi == 0.0
    assert out.dF_deta == 0.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n": 0},
        {"X": -1.0},
        {"eta": 1.1},
        {"C": 0.0},
        {"phi": (lambda x, e: 1.0,)},
    ],
)
def test_eq_5_15_bad_inputs_fail_closed(kwargs):
    args = dict(
        n=1,
        X=0.5,
        eta=0.0,
        U_n=lambda x, e: 1.0,
        dU_n_deta=lambda x, e: 0.0,
        phi=(lambda x, e: 1.0, lambda x, e: 0.5),
        omega_prev_over_x=lambda x, e: 0.0,
        h=0.005,
        C=2.0,
    )
    args.update(kwargs)
    with pytest.raises((ValueError, TypeError)):
        reconstruct_eq_5_15(**args)


def test_eq_5_15_nonfinite_profile_fails_closed():
    with pytest.raises(ValueError):
        reconstruct_eq_5_15(
            1,
            0.5,
            0.0,
            lambda x, e: math.nan,
            lambda x, e: 0.0,
            (lambda x, e: 1.0, lambda x, e: 0.5),
            lambda x, e: 0.0,
            h=0.005,
            C=2.0,
        )
