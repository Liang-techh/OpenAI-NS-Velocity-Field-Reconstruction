import math

import pytest

from openai_ns_reconstruction.background_lower_history_source import (
    preceding_diffusion_from_second_jet,
)
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
    preceding_diffusion_parameter_jet,
)


H = 0.005

# F(X,eta) = sum a_ij X^i eta^j.  Degrees (2,3) exercise exactly the
# third-mixed derivatives needed by the analytic eta derivative.
_COEFFICIENTS = {
    (0, 0): 0.7,
    (0, 1): -0.2,
    (0, 2): 0.11,
    (0, 3): -0.04,
    (1, 0): 0.3,
    (1, 1): 0.17,
    (1, 2): -0.09,
    (1, 3): 0.025,
    (2, 0): -0.08,
    (2, 1): 0.06,
    (2, 2): 0.031,
    (2, 3): -0.012,
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
    return ProfileThirdMixedJet(
        value=_derivative(X, eta, 0, 0),
        radial=_derivative(X, eta, 1, 0),
        radial2=_derivative(X, eta, 2, 0),
        parameter=_derivative(X, eta, 0, 1),
        radial_parameter=_derivative(X, eta, 1, 1),
        parameter2=_derivative(X, eta, 0, 2),
        radial2_parameter=_derivative(X, eta, 2, 1),
        radial_parameter2=_derivative(X, eta, 1, 2),
        parameter3=_derivative(X, eta, 0, 3),
    )


@pytest.mark.parametrize("order", [1, 2, 4])
@pytest.mark.parametrize("X", [0.0, 0.23, 1.1])
@pytest.mark.parametrize("power", [-1.005, -0.505])
def test_parameter_jet_matches_independent_value_path_finite_difference(
    order, X, power
):
    eta = 0.31
    result = preceding_diffusion_parameter_jet(
        H, power, order, X, eta, _jet(X, eta)
    )

    expected_value = preceding_diffusion_from_second_jet(
        H, power, order, X, eta, _jet(X, eta).second()
    )
    assert result.value == pytest.approx(expected_value, rel=0.0, abs=1e-13)

    epsilon = 2.0e-6
    plus = preceding_diffusion_from_second_jet(
        H, power, order, X, eta + epsilon, _jet(X, eta + epsilon).second()
    )
    minus = preceding_diffusion_from_second_jet(
        H, power, order, X, eta - epsilon, _jet(X, eta - epsilon).second()
    )
    finite_difference = (plus - minus) / (2.0 * epsilon)
    assert result.parameter == pytest.approx(
        finite_difference, rel=2.0e-7, abs=2.0e-8
    )


def test_all_three_new_third_mixed_entries_are_used():
    X = 0.73
    eta = 0.28
    base = _jet(X, eta)
    reference = preceding_diffusion_parameter_jet(
        H, -1.005, 3, X, eta, base
    ).parameter

    names = ("radial2_parameter", "radial_parameter2", "parameter3")
    for name in names:
        values = base.__dict__.copy()
        values[name] += 0.37
        changed = preceding_diffusion_parameter_jet(
            H, -1.005, 3, X, eta, ProfileThirdMixedJet(**values)
        ).parameter
        assert abs(changed - reference) > 1.0e-6


def test_parameter_jet_fails_closed_on_missing_stronger_jet():
    with pytest.raises(TypeError, match="ProfileThirdMixedJet"):
        preceding_diffusion_parameter_jet(H, -1.005, 1, 0.2, 0.1, object())

    with pytest.raises(ValueError, match="positive integer"):
        preceding_diffusion_parameter_jet(H, -1.005, 0, 0.2, 0.1, _jet(0.2, 0.1))

    with pytest.raises(ValueError, match="nonnegative"):
        preceding_diffusion_parameter_jet(H, -1.005, 1, -0.1, 0.1, _jet(0.2, 0.1))
