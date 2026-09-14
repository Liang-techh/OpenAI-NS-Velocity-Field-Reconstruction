import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_second_parameter_jet import (
    preceding_diffusion_second_parameter_jet,
)
from openai_ns_reconstruction.background_preceding_diffusion_third_parameter_jet import (
    hierarchy_owned_angular_preceding_diffusion_third_parameter_jet,
    preceding_diffusion_third_parameter_jet,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_phi_fifth_mixed import (
    Section5LowerHistoryPhiFifthMixedHierarchy,
    Section5PhiFifthMixedCoefficientJetSource,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


H = 0.005
C = 2.75

# Degree (2,5) makes all three newly required total-order-five entries nonzero.
_COEFFICIENTS = {
    (0, 0): 0.7,
    (0, 1): -0.2,
    (0, 2): 0.11,
    (0, 3): -0.04,
    (0, 4): 0.013,
    (0, 5): -0.003,
    (1, 0): 0.3,
    (1, 1): 0.17,
    (1, 2): -0.09,
    (1, 3): 0.025,
    (1, 4): -0.008,
    (1, 5): 0.0017,
    (2, 0): -0.08,
    (2, 1): 0.06,
    (2, 2): 0.031,
    (2, 3): -0.012,
    (2, 4): 0.004,
    (2, 5): -0.0009,
}


def _falling(power, derivative):
    out = 1
    for j in range(derivative):
        out *= power - j
    return out


def _derivative(X, eta, radial, parameter):
    return math.fsum(
        coefficient
        * _falling(i, radial)
        * _falling(j, parameter)
        * X ** (i - radial)
        * eta ** (j - parameter)
        for (i, j), coefficient in _COEFFICIENTS.items()
        if i >= radial and j >= parameter
    )


def _jet(X, eta):
    return ProfileFifthMixedJet(
        value=_derivative(X, eta, 0, 0),
        radial=_derivative(X, eta, 1, 0),
        radial2=_derivative(X, eta, 2, 0),
        parameter=_derivative(X, eta, 0, 1),
        radial_parameter=_derivative(X, eta, 1, 1),
        parameter2=_derivative(X, eta, 0, 2),
        radial2_parameter=_derivative(X, eta, 2, 1),
        radial_parameter2=_derivative(X, eta, 1, 2),
        parameter3=_derivative(X, eta, 0, 3),
        radial2_parameter2=_derivative(X, eta, 2, 2),
        radial_parameter3=_derivative(X, eta, 1, 3),
        parameter4=_derivative(X, eta, 0, 4),
        radial2_parameter3=_derivative(X, eta, 2, 3),
        radial_parameter4=_derivative(X, eta, 1, 4),
        parameter5=_derivative(X, eta, 0, 5),
    )


@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("X", [0.0, 0.23, 1.1])
@pytest.mark.parametrize("power", [-1.005, -0.505])
def test_third_parameter_jet_preserves_lower_rows_and_matches_fd_oracle(
    order, X, power
):
    eta = 0.31
    result = preceding_diffusion_third_parameter_jet(
        H, power, order, X, eta, _jet(X, eta)
    )
    lower = preceding_diffusion_second_parameter_jet(
        H, power, order, X, eta, _jet(X, eta).fourth()
    )
    assert result.value == lower.value
    assert result.parameter == lower.parameter
    assert result.parameter2 == lower.parameter2

    epsilon = 2.0e-5
    plus = preceding_diffusion_second_parameter_jet(
        H, power, order, X, eta + epsilon, _jet(X, eta + epsilon).fourth()
    ).parameter2
    minus = preceding_diffusion_second_parameter_jet(
        H, power, order, X, eta - epsilon, _jet(X, eta - epsilon).fourth()
    ).parameter2
    finite_difference = (plus - minus) / (2.0 * epsilon)
    assert result.parameter3 == pytest.approx(
        finite_difference, rel=1.2e-6, abs=1.2e-7
    )


def test_all_three_new_fifth_mixed_entries_are_used():
    X = 0.73
    eta = 0.28
    base = _jet(X, eta)
    reference = preceding_diffusion_third_parameter_jet(
        H, -1.005, 3, X, eta, base
    ).parameter3

    names = ("radial2_parameter3", "radial_parameter4", "parameter5")
    for name in names:
        values = base.__dict__.copy()
        values[name] += 0.37
        changed = preceding_diffusion_third_parameter_jet(
            H, -1.005, 3, X, eta, ProfileFifthMixedJet(**values)
        ).parameter3
        assert abs(changed - reference) > 1.0e-6


def test_third_parameter_jet_rejects_weaker_and_nonfinite_new_entries():
    with pytest.raises(TypeError, match="ProfileFifthMixedJet"):
        preceding_diffusion_third_parameter_jet(
            H, -1.005, 1, 0.2, 0.1, _jet(0.2, 0.1).fourth()
        )

    values = _jet(0.2, 0.1).__dict__.copy()
    values["parameter5"] = math.nan
    with pytest.raises(ValueError, match="parameter5 must be finite"):
        preceding_diffusion_third_parameter_jet(
            H, -1.005, 1, 0.2, 0.1, ProfileFifthMixedJet(**values)
        )


_MOMENT_POLYS = np.array(
    [
        [0.20, 0.10, 0.03, 0.004, 0.0007, 0.00011],
        [-0.15, 0.00, 0.02, -0.003, 0.0005, -0.00009],
        [0.11, -0.03, 0.00, 0.002, -0.0004, 0.00008],
        [0.25, 0.04, 0.01, -0.005, 0.0006, -0.00007],
        [-0.12, 0.05, 0.006, 0.001, -0.0003, 0.00006],
    ],
    dtype=float,
)


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
            for derivative in range(6)
        ]
    )


def _patch_jets(_eta):
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0, 0.0])


def _functional_repair():
    repair = Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )
    return Lemma52RepairedProfileAdapter(
        repair,
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )


def _base_phi_fifth_mixed_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    coefficients = [1.0, 0.04, 0.05, 0.01, 0.002, 0.0003]
    rows = [
        e * _poly_derivative(coefficients, eta, derivative)
        for derivative in range(6)
    ]
    value, parameter, parameter2, parameter3, parameter4, parameter5 = rows
    return ProfileFifthMixedJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=parameter2,
        radial2_parameter=parameter,
        radial_parameter2=-parameter2,
        parameter3=parameter3,
        radial2_parameter2=parameter2,
        radial_parameter3=-parameter3,
        parameter4=parameter4,
        radial2_parameter3=parameter3,
        radial_parameter4=-parameter4,
        parameter5=parameter5,
    )


def _base_u_fifth_mixed_jet(X, eta):
    e = math.exp(-float(X))
    eta = float(eta)
    return AxialFifthMixedJet(
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
        radial2_parameter3=0.0,
        radial_parameter4=0.0,
        parameter5=0.0,
    )


def _leading_source():
    phi = _base_phi_fifth_mixed_jet
    axial = _base_u_fifth_mixed_jet
    return Section5PhiFifthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: phi(X, eta).fourth().third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: axial(X, eta).fourth().third(),
        provenance="analytic leading fifth-mixed test fixture; Issue #1 remains upstream",
        axial_fourth_mixed_jet_provider=lambda X, eta: axial(X, eta).fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: phi(X, eta).fourth().third(),
        axial_fifth_mixed_jet_provider=axial,
        phi_fourth_mixed_jet_provider=lambda X, eta: phi(X, eta).fourth(),
        phi_fifth_mixed_jet_provider=phi,
    )


def _hierarchy():
    repaired = Section5PhiFifthMixedCoefficientJetSource.from_lemma52_repair_full_strong_mixed(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_fifth_mixed_jet=_base_phi_fifth_mixed_jet,
        base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
    )
    return Section5LowerHistoryPhiFifthMixedHierarchy(
        H,
        C,
        (_leading_source(), repaired),
        quadrature_points=64,
    )


def test_hierarchy_owned_angular_bridge_consumes_real_compact_repaired_phi_fifth_jet():
    hierarchy = _hierarchy()
    X = 0.5 * 4.17**2  # inside the first active E-repair bump
    eta = -0.23

    repaired_jet = hierarchy.phi_fifth_mixed_jet(1, X, eta)
    assert repaired_jet != _base_phi_fifth_mixed_jet(X, eta)

    expected = preceding_diffusion_third_parameter_jet(
        H,
        -1.0 - H,
        2,
        X,
        eta,
        repaired_jet,
    )
    actual = hierarchy_owned_angular_preceding_diffusion_third_parameter_jet(
        hierarchy,
        2,
        X,
        eta,
    )
    assert actual == expected


def test_hierarchy_owned_bridge_rejects_missing_strict_lower_order():
    hierarchy = Section5LowerHistoryPhiFifthMixedHierarchy(
        H,
        C,
        (_leading_source(),),
    )
    with pytest.raises(ValueError, match="missing one or more strict-lower"):
        hierarchy_owned_angular_preceding_diffusion_third_parameter_jet(
            hierarchy,
            2,
            0.0,
            0.1,
        )
