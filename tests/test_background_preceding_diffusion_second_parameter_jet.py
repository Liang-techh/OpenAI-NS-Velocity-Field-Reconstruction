import math

import pytest

from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    preceding_diffusion_parameter_jet,
)
from openai_ns_reconstruction.background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
    preceding_diffusion_second_parameter_jet,
)


H = 0.005

# Degree (2,4) makes every newly required fourth mixed derivative nontrivial.
_COEFFICIENTS = {
    (0, 0): 0.7,
    (0, 1): -0.2,
    (0, 2): 0.11,
    (0, 3): -0.04,
    (0, 4): 0.013,
    (1, 0): 0.3,
    (1, 1): 0.17,
    (1, 2): -0.09,
    (1, 3): 0.025,
    (1, 4): -0.008,
    (2, 0): -0.08,
    (2, 1): 0.06,
    (2, 2): 0.031,
    (2, 3): -0.012,
    (2, 4): 0.004,
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
    return ProfileFourthMixedJet(
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
    )


@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("X", [0.0, 0.23, 1.1])
@pytest.mark.parametrize("power", [-1.005, -0.505])
def test_second_parameter_jet_matches_landed_first_jet_and_fd_oracle(
    order, X, power
):
    eta = 0.31
    result = preceding_diffusion_second_parameter_jet(
        H, power, order, X, eta, _jet(X, eta)
    )
    first = preceding_diffusion_parameter_jet(
        H, power, order, X, eta, _jet(X, eta).third()
    )
    assert result.value == first.value
    assert result.parameter == first.parameter

    epsilon = 2.0e-5
    plus = preceding_diffusion_parameter_jet(
        H, power, order, X, eta + epsilon, _jet(X, eta + epsilon).third()
    ).parameter
    minus = preceding_diffusion_parameter_jet(
        H, power, order, X, eta - epsilon, _jet(X, eta - epsilon).third()
    ).parameter
    finite_difference = (plus - minus) / (2.0 * epsilon)
    assert result.parameter2 == pytest.approx(
        finite_difference, rel=8.0e-7, abs=8.0e-8
    )


def test_all_three_new_fourth_mixed_entries_are_used():
    X = 0.73
    eta = 0.28
    base = _jet(X, eta)
    reference = preceding_diffusion_second_parameter_jet(
        H, -1.005, 3, X, eta, base
    ).parameter2

    names = ("radial2_parameter2", "radial_parameter3", "parameter4")
    for name in names:
        values = base.__dict__.copy()
        values[name] += 0.37
        changed = preceding_diffusion_second_parameter_jet(
            H, -1.005, 3, X, eta, ProfileFourthMixedJet(**values)
        ).parameter2
        assert abs(changed - reference) > 1.0e-6


def test_second_parameter_jet_fails_closed_on_weaker_or_invalid_input():
    with pytest.raises(TypeError, match="ProfileFourthMixedJet"):
        preceding_diffusion_second_parameter_jet(
            H, -1.005, 1, 0.2, 0.1, _jet(0.2, 0.1).third()
        )

    values = _jet(0.2, 0.1).__dict__.copy()
    values["parameter4"] = math.nan
    with pytest.raises(ValueError, match="parameter4 must be finite"):
        ProfileFourthMixedJet(**values)
