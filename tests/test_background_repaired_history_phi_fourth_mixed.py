import math
from dataclasses import replace

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_fifth_mixed_jets import (
    Lemma52RepairedFifthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_phi_fourth_mixed_jets import (
    Lemma52RepairedPhiFourthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
    preceding_diffusion_second_parameter_jet,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_fifth_mixed import (
    Section5FifthMixedCoefficientJetSource,
)
from openai_ns_reconstruction.background_repaired_history_phi_fourth_mixed import (
    Section5LowerHistoryPhiFourthMixedHierarchy,
    Section5PhiFourthMixedCoefficientJetSource,
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
        [0.20, 0.10, 0.03, 0.004, 0.0007, 0.00011],
        [-0.15, 0.00, 0.02, -0.003, 0.0005, -0.00009],
        [0.11, -0.03, 0.00, 0.002, -0.0004, 0.00008],
        [0.25, 0.04, 0.01, -0.005, 0.0006, -0.00007],
        [-0.12, 0.05, 0.006, 0.001, -0.0003, 0.00006],
    ],
    dtype=float,
)


def _falling(power, derivative):
    if derivative > power:
        return 0
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
            for derivative in range(6)
        ]
    )


def _patch_jets(_eta):
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0, 0.0])


def _functional_repair():
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )


def _base_phi_fourth_mixed_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    coeff = [1.0, 0.04, 0.05, 0.01, 0.002]
    rows = [e * _poly_derivative(coeff, eta, derivative) for derivative in range(5)]
    value, parameter, parameter2, parameter3, parameter4 = rows
    return ProfileFourthMixedJet(
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


def _leading_strong_source():
    return Section5PhiFourthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance=(
            "analytic phi-fourth/U-fifth leading test fixture; Issue #1 remains upstream"
        ),
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third(),
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
        phi_fourth_mixed_jet_provider=_base_phi_fourth_mixed_jet,
    )


def _hierarchy():
    repaired = (
        Section5PhiFourthMixedCoefficientJetSource.from_lemma52_repair_full_strong_mixed(
            1,
            C=C,
            profile=_functional_repair(),
            base_phi_fourth_mixed_jet=_base_phi_fourth_mixed_jet,
            base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
        )
    )
    return Section5LowerHistoryPhiFourthMixedHierarchy(
        H,
        C,
        (_leading_strong_source(), repaired),
        quadrature_points=64,
    )


def test_hierarchy_owns_repaired_phi_fourth_and_exact_lower_projections():
    functional = _functional_repair()
    hierarchy = _hierarchy()
    X = 0.5 * 4.17**2
    eta = -0.23

    independent = Lemma52RepairedPhiFourthMixedJetAdapter(
        functional,
        C,
        _base_phi_fourth_mixed_jet,
    )
    expected = independent.phi_fourth_mixed_jet(X, eta)
    actual = hierarchy.phi_fourth_mixed_jet(1, X, eta)

    assert actual == expected
    assert actual.third() == hierarchy.phi_third_mixed_jet(1, X, eta)
    assert actual.third().second() == hierarchy.phi_second_jet(1, X, eta)


def test_same_strong_hierarchy_preserves_repaired_fifth_mixed_u_ownership():
    functional = _functional_repair()
    hierarchy = _hierarchy()
    X = 0.5 * 2.17**2
    eta = 0.19

    independent = Lemma52RepairedFifthMixedJetAdapter(
        functional,
        _base_u_fifth_mixed_jet,
    )
    assert hierarchy.axial_fifth_mixed_jet(1, X, eta) == independent.U_fifth_mixed_jet(
        X, eta
    )


def test_owned_phi_fourth_drives_existing_second_eta_preceding_diffusion():
    hierarchy = _hierarchy()
    X = 0.5 * 4.17**2
    eta = 0.21
    power = -1.0 - H

    actual = preceding_diffusion_second_parameter_jet(
        H,
        power,
        2,
        X,
        eta,
        hierarchy.phi_fourth_mixed_jet(1, X, eta),
    )
    expected = preceding_diffusion_second_parameter_jet(
        H,
        power,
        2,
        X,
        eta,
        Lemma52RepairedPhiFourthMixedJetAdapter(
            _functional_repair(), C, _base_phi_fourth_mixed_jet
        ).phi_fourth_mixed_jet(X, eta),
    )
    assert actual == expected


def test_phi_fourth_query_fails_closed_for_previous_fifth_mixed_source():
    weak = Section5FifthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance="previous strong source without fourth-mixed phi ownership",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third(),
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
    )
    hierarchy = Section5LowerHistoryPhiFourthMixedHierarchy(H, C, (weak,))

    with pytest.raises(ValueError, match="fourth-mixed phi"):
        hierarchy.phi_fourth_mixed_jet(0, 0.0, 0.1)


def test_strong_hierarchy_rejects_incoherent_fourth_to_third_projection():
    def bad_third(X, eta):
        third = _base_phi_fourth_mixed_jet(X, eta).third()
        return replace(third, parameter3=third.parameter3 + 1.0)

    bad = Section5PhiFourthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance="deliberately incoherent phi fourth/third test fixture",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=bad_third,
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
        phi_fourth_mixed_jet_provider=_base_phi_fourth_mixed_jet,
    )
    hierarchy = Section5LowerHistoryPhiFourthMixedHierarchy(H, C, (bad,))

    with pytest.raises(ValueError, match="incoherent"):
        hierarchy.phi_third_mixed_jet(0, 0.5 * 4.17**2, 0.2)


def test_lemma52_full_strong_constructor_rejects_order_zero():
    with pytest.raises(ValueError, match="positive integer"):
        Section5PhiFourthMixedCoefficientJetSource.from_lemma52_repair_full_strong_mixed(
            0,
            C=C,
            profile=_functional_repair(),
            base_phi_fourth_mixed_jet=_base_phi_fourth_mixed_jet,
            base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
        )
