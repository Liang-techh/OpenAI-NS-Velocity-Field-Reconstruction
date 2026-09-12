import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_fourth_mixed_jets import (
    Lemma52RepairedFourthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_regular_flux_second_jets import (
    AxialThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history import (
    Section5CoefficientJetSource,
    Section5LowerHistoryJetHierarchy,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


C = 2.75
H = 0.005


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


_MOMENT_POLYS = np.array(
    [
        [0.20, 0.10, 0.03, 0.004, 0.0007],
        [-0.15, 0.00, 0.02, -0.003, 0.0005],
        [0.11, -0.03, 0.00, 0.002, -0.0004],
        [0.25, 0.04, 0.01, -0.005, 0.0006],
        [-0.12, 0.05, 0.006, 0.001, -0.0003],
    ],
    dtype=float,
)


def _falling(power, derivative):
    out = 1
    for j in range(derivative):
        out *= power - j
    return out


def _poly_derivative(coefficients, eta, derivative):
    return math.fsum(
        float(coefficient)
        * _falling(power, derivative)
        * float(eta) ** (power - derivative)
        for power, coefficient in enumerate(coefficients)
        if power >= derivative
    )


def _moment_jets(eta):
    return np.vstack(
        [
            np.array(
                [_poly_derivative(row, eta, derivative) for row in _MOMENT_POLYS],
                dtype=float,
            )
            for derivative in range(5)
        ]
    )


def _patch_jets(_eta):
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0])


def _functional_repair():
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )


def _base_phi_second_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    angular = 1.0 + 0.1 * eta * eta
    value = C * e * angular
    parameter = C * e * 0.2 * eta
    return ProfileSecondJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=C * e * 0.2,
    )


def _base_u_fourth_mixed_jet(X, eta):
    e = math.exp(-float(X))
    eta = float(eta)
    return AxialFourthMixedJet(
        value=eta * e,
        radial=-eta * e,
        radial2=eta * e,
        parameter=e,
        radial_parameter=-e,
        parameter2=0.0,
        radial2_parameter=e,
        radial_parameter2=0.0,
        parameter3=0.0,
        radial2_parameter2=0.0,
        radial_parameter3=0.0,
        parameter4=0.0,
    )


def _leading_third_mixed_jet(X, eta):
    return _base_u_fourth_mixed_jet(X, eta).third()


def _hierarchy_with_fourth():
    leading = Section5CoefficientJetSource(
        0,
        _base_phi_second_jet,
        _leading_third_mixed_jet,
        "analytic leading-order test fixture; Issue #1 remains upstream",
    )
    repaired = Section5CoefficientJetSource.from_lemma52_repair_fourth_mixed(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_second_jet=_base_phi_second_jet,
        base_U_fourth_mixed_jet=_base_u_fourth_mixed_jet,
    )
    return Section5LowerHistoryJetHierarchy(
        H,
        C,
        (leading, repaired),
        quadrature_points=64,
    )


def test_hierarchy_owns_actual_repaired_fourth_mixed_u_and_derives_beta_third_mixed():
    hierarchy = _hierarchy_with_fourth()
    functional = _functional_repair()
    X = 0.5 * 2.17**2
    eta = 0.23

    independent = Lemma52RepairedFourthMixedJetAdapter(
        functional, _base_u_fourth_mixed_jet
    )
    expected_u = independent.U_fourth_mixed_jet(X, eta)
    actual_u = hierarchy.axial_fourth_mixed_jet(1, X, eta)

    assert actual_u == expected_u
    assert hierarchy.axial_third_mixed_jet(1, X, eta) == expected_u.third()

    expected_beta = independent.beta_third_mixed_jet(
        H, 1, X, eta, quadrature_points=64
    )
    actual_beta = hierarchy.beta_third_mixed_jet(1, X, eta)
    assert actual_beta == expected_beta

    # The ordinary lower-history beta path is the exact projection of the same
    # hierarchy-owned fourth-mixed U source, not a second independent record.
    beta_second = hierarchy.beta_second_jet(1, X, eta)
    assert beta_second.value == actual_beta.value
    assert beta_second.radial2 == actual_beta.radial2
    assert beta_second.parameter2 == actual_beta.parameter2


def test_fourth_mixed_path_fails_closed_when_coefficient_only_owns_old_third_layer():
    leading = Section5CoefficientJetSource(
        0,
        _base_phi_second_jet,
        _leading_third_mixed_jet,
        "leading fixture",
    )
    repaired_old = Section5CoefficientJetSource.from_lemma52_repair(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_second_jet=_base_phi_second_jet,
        base_U_third_mixed_jet=_leading_third_mixed_jet,
    )
    hierarchy = Section5LowerHistoryJetHierarchy(H, C, (leading, repaired_old))

    with pytest.raises(ValueError, match="fourth-mixed U"):
        hierarchy.axial_fourth_mixed_jet(1, 0.5 * 2.17**2, 0.1)
    with pytest.raises(ValueError, match="fourth-mixed U"):
        hierarchy.beta_third_mixed_jet(1, 0.5 * 2.17**2, 0.1)


def test_fourth_mixed_constructor_rejects_order_zero_and_non_callable_base_provider():
    with pytest.raises(ValueError):
        Section5CoefficientJetSource.from_lemma52_repair_fourth_mixed(
            0,
            C=C,
            profile=_functional_repair(),
            base_phi_second_jet=_base_phi_second_jet,
            base_U_fourth_mixed_jet=_base_u_fourth_mixed_jet,
        )
    with pytest.raises(TypeError, match="base_U_fourth_mixed_jet"):
        Section5CoefficientJetSource.from_lemma52_repair_fourth_mixed(
            1,
            C=C,
            profile=_functional_repair(),
            base_phi_second_jet=_base_phi_second_jet,
            base_U_fourth_mixed_jet=None,
        )
