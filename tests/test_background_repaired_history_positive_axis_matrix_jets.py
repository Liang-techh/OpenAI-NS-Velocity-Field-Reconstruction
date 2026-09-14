import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_positive_axis_base import (
    hierarchy_owned_positive_axis_matrices,
)
from openai_ns_reconstruction.background_repaired_history_positive_axis_matrix_jets import (
    hierarchy_owned_positive_axis_matrix_second_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
    Section5SixthMixedCoefficientJetSource,
)


H = 0.005
C = 2.75


def _poly_derivatives(coefficients, eta, max_order):
    values = []
    for derivative in range(max_order + 1):
        total = 0.0
        for power in range(derivative, len(coefficients)):
            factor = math.factorial(power) / math.factorial(power - derivative)
            total += coefficients[power] * factor * eta ** (power - derivative)
        values.append(total)
    return values


def _phi(X, eta):
    g = _poly_derivatives((1.0, 0.04, 0.015, -0.006, 0.002, -0.0005), eta, 5)
    radial0 = 0.31 * (1.0 + 0.17 * X + 0.03 * X * X)
    radial1 = 0.31 * (0.17 + 0.06 * X)
    radial2 = 0.31 * 0.06
    return ProfileFifthMixedJet(
        value=radial0 * g[0],
        radial=radial1 * g[0],
        radial2=radial2 * g[0],
        parameter=radial0 * g[1],
        radial_parameter=radial1 * g[1],
        parameter2=radial0 * g[2],
        radial2_parameter=radial2 * g[1],
        radial_parameter2=radial1 * g[2],
        parameter3=radial0 * g[3],
        radial2_parameter2=radial2 * g[2],
        radial_parameter3=radial1 * g[3],
        parameter4=radial0 * g[4],
        radial2_parameter3=radial2 * g[3],
        radial_parameter4=radial1 * g[4],
        parameter5=radial0 * g[5],
    )


def _axial(X, eta):
    g = _poly_derivatives(
        (1.0, -0.03, 0.01, 0.004, -0.0015, 0.0004, -0.00008),
        eta,
        6,
    )
    radial0 = 0.22 * (1.0 + 0.11 * X + 0.025 * X * X)
    radial1 = 0.22 * (0.11 + 0.05 * X)
    radial2 = 0.22 * 0.05
    return AxialSixthMixedJet(
        value=radial0 * g[0],
        radial=radial1 * g[0],
        radial2=radial2 * g[0],
        parameter=radial0 * g[1],
        radial_parameter=radial1 * g[1],
        parameter2=radial0 * g[2],
        radial2_parameter=radial2 * g[1],
        radial_parameter2=radial1 * g[2],
        parameter3=radial0 * g[3],
        radial2_parameter2=radial2 * g[2],
        radial_parameter3=radial1 * g[3],
        parameter4=radial0 * g[4],
        radial2_parameter3=radial2 * g[3],
        radial_parameter4=radial1 * g[4],
        parameter5=radial0 * g[5],
        radial2_parameter4=radial2 * g[4],
        radial_parameter5=radial1 * g[5],
        parameter6=radial0 * g[6],
    )


def _source0():
    return Section5SixthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _phi(X, eta).fourth().third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _axial(X, eta).fifth().fourth().third(),
        provenance="analytic PositiveAxis matrix-jet regression fixture; not paper coefficient data",
        normalization_C=None,
        axial_fourth_mixed_jet_provider=lambda X, eta: _axial(X, eta).fifth().fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _phi(X, eta).fourth().third(),
        axial_fifth_mixed_jet_provider=lambda X, eta: _axial(X, eta).fifth(),
        phi_fourth_mixed_jet_provider=lambda X, eta: _phi(X, eta).fourth(),
        phi_fifth_mixed_jet_provider=_phi,
        axial_sixth_mixed_jet_provider=_axial,
    )


def _hierarchy():
    return Section5LowerHistorySixthMixedHierarchy(
        H,
        C,
        (_source0(),),
        quadrature_points=8,
    )


class _Series2:
    """Test-only normalized eta Taylor coefficients through degree two."""

    def __init__(self, c0, c1=0.0, c2=0.0):
        self.c = (float(c0), float(c1), float(c2))

    @classmethod
    def derivatives(cls, value, d1, d2):
        return cls(value, d1, 0.5 * d2)

    @classmethod
    def coerce(cls, value):
        return value if isinstance(value, cls) else cls(value)

    def __add__(self, other):
        other = self.coerce(other)
        return _Series2(*(a + b for a, b in zip(self.c, other.c)))

    __radd__ = __add__

    def __neg__(self):
        return _Series2(*(-a for a in self.c))

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        a, b = self.c, other.c
        return _Series2(
            a[0] * b[0],
            a[0] * b[1] + a[1] * b[0],
            a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        a, b = self.c, other.c
        q0 = a[0] / b[0]
        q1 = (a[1] - b[1] * q0) / b[0]
        q2 = (a[2] - b[1] * q1 - b[2] * q0) / b[0]
        return _Series2(q0, q1, q2)

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    @property
    def d1(self):
        return self.c[1]

    @property
    def d2(self):
        return 2.0 * self.c[2]


def _matrix_series_oracle(hierarchy, order, xi, eta):
    """Independent normalized-series replay of the displayed A0/A1 formulas."""

    X = xi * xi
    phi = hierarchy.phi_fifth_mixed_jet(0, X, eta)
    axial = hierarchy.axial_sixth_mixed_jet(0, X, eta)
    beta = hierarchy.beta_fifth_mixed_jet(0, X, eta)

    p = _Series2.derivatives(phi.value, phi.parameter, phi.parameter2)
    px = _Series2.derivatives(
        phi.radial,
        phi.radial_parameter,
        phi.radial_parameter2,
    )
    pe = _Series2.derivatives(phi.parameter, phi.parameter2, phi.parameter3)
    u = _Series2.derivatives(axial.value, axial.parameter, axial.parameter2)
    ux = _Series2.derivatives(
        axial.radial,
        axial.radial_parameter,
        axial.radial_parameter2,
    )
    ue = _Series2.derivatives(axial.parameter, axial.parameter2, axial.parameter3)
    beta_j = _Series2.derivatives(beta.value, beta.parameter, beta.parameter2)

    e = _Series2.derivatives(eta, 1.0, 0.0)
    h = hierarchy.h
    a = 0.5 + h
    d = 0.5 - h
    lam = 2.0 * order * h
    edge = 1.0 - e * e
    ell = 1.0 - 2.0 * h * e * e
    inv_c2 = 1.0 / (hierarchy.C * hierarchy.C)
    M = 1.0 - 2.0 * e * u
    R = M / ell + beta_j
    Gphi = X * px + p
    Gu = X * ux

    def transported(power, value, parameter, radial):
        return (
            2.0 * e * power * value
            + edge * parameter
            - 2.0 * e * X * radial
        ) / ell

    A0 = {}
    A0[3, 0] = 4.0 * xi * inv_c2 * p
    A0[4, 0] = 2.0 * (beta_j - (-a - 0.5 + lam) * M / ell)
    A0[4, 1] = 2.0 * (
        transported(-a - 0.5, p, pe, px)
        + 2.0 * e * (a - lam) * Gphi / ell
    )
    A0[4, 2] = -4.0 * e * (d + lam) * Gphi / ell
    A0[4, 4] = xi * R
    A0[5, 0] = -8.0 * e * X * inv_c2 * p / ell
    A0[5, 1] = 2.0 * (
        -(-a + lam) * M / ell
        + transported(-a, u, ue, ux)
        + 2.0 * e * (a - lam) * Gu / ell
    )
    A0[5, 2] = -4.0 * e * (d + lam) * Gu / ell
    A0[5, 3] = 4.0 * e * (-2.0 * a + lam) / ell
    A0[5, 5] = xi * R

    H = d * e + edge * u
    A1 = {}
    A1[4, 0] = 2.0 * H / ell
    A1[4, 1] = -2.0 * edge * Gphi / ell
    A1[4, 2] = -2.0 * edge * Gphi / ell
    A1[5, 1] = 2.0 * (H - edge * Gu) / ell
    A1[5, 2] = -2.0 * edge * Gu / ell
    A1[5, 3] = 2.0 * edge / ell

    out = []
    for entries in (A0, A1):
        d1 = np.zeros((6, 6), dtype=float)
        d2 = np.zeros((6, 6), dtype=float)
        for index, series in entries.items():
            d1[index] = series.d1
            d2[index] = series.d2
        out.extend((d1, d2))
    return tuple(out)


def test_matrix_eta2_value_rows_delegate_exactly_to_landed_hierarchy_values():
    hierarchy = _hierarchy()
    order = 1
    xi = 0.37
    eta = -0.21
    actual = hierarchy_owned_positive_axis_matrix_second_parameter_jet(
        hierarchy,
        order,
        xi,
        eta,
    )
    landed = hierarchy_owned_positive_axis_matrices(hierarchy, order, xi, eta)
    np.testing.assert_array_equal(actual.A0, landed.A0)
    np.testing.assert_array_equal(actual.A1, landed.A1)


def test_matrix_eta2_matches_independent_normalized_taylor_series_oracle():
    hierarchy = _hierarchy()
    order = 2
    xi = 0.41
    eta = 0.26
    actual = hierarchy_owned_positive_axis_matrix_second_parameter_jet(
        hierarchy,
        order,
        xi,
        eta,
    )
    A0_d1, A0_d2, A1_d1, A1_d2 = _matrix_series_oracle(
        hierarchy,
        order,
        xi,
        eta,
    )
    np.testing.assert_allclose(actual.A0_parameter, A0_d1, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(actual.A0_parameter2, A0_d2, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(actual.A1_parameter, A1_d1, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(actual.A1_parameter2, A1_d2, rtol=2e-13, atol=2e-13)


def test_matrix_eta2_is_axis_regular_read_only_and_finite():
    actual = hierarchy_owned_positive_axis_matrix_second_parameter_jet(
        _hierarchy(),
        1,
        0.0,
        0.17,
    )
    for matrix in (
        actual.A0,
        actual.A0_parameter,
        actual.A0_parameter2,
        actual.A1,
        actual.A1_parameter,
        actual.A1_parameter2,
    ):
        assert np.all(np.isfinite(matrix))
        assert not matrix.flags.writeable


def test_matrix_eta2_fails_closed_on_wrong_hierarchy_type():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_positive_axis_matrix_second_parameter_jet(
            object(),
            1,
            0.0,
            0.0,
        )
