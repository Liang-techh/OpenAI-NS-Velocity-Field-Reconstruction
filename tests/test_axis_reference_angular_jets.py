import math

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.axis_reference_angular_jets import (
    actual_schedule_angular_reference_jets,
    axis_weight,
)
from openai_ns_reconstruction.natural_axis import D, H, chi
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def jets():
    return actual_schedule_angular_reference_jets(_schedule_data(), 0.05)


def test_angular_jets_use_actual_schedule_sigma_and_epsilon(jets) -> None:
    assert jets.reference.sigma == jets.analytic_inputs.neighborhood.sigma
    assert jets.epsilon == jets.analytic_inputs.neighborhood.radius / 2.0
    assert jets.epsilon > 0.0
    assert jets.paper_exact is False


def test_zeroth_parameter_jet_matches_landed_reference_pair(jets) -> None:
    for eta in (-1.1, -0.37, 0.0, 0.42, 1.1):
        for n in range(7):
            assert jets.parameter_jet(n, 0, eta) == pytest.approx(
                jets.reference.phi_coefficient(n, eta), rel=3e-14, abs=1e-300
            )


def test_first_parameter_jet_matches_independent_closed_form(jets) -> None:
    eta = 0.31
    h = jets.reference.data.h
    j = jets.reference.j
    sigma = jets.reference.sigma
    hh = H(h, j, eta)
    hp = D(h) + 4.0 - 12.0 * eta * eta - 2.0 * j * eta
    q = hh * hh
    chi_value = chi(h, j, sigma, eta)
    chi_prime = (2.0 * hh * hp * sigma * sigma) / (q + sigma * sigma) ** 2
    c = -0.5 * chi_value
    c_prime = -0.5 * chi_prime

    assert jets.parameter_jet(0, 1, eta) == 0.0
    for n in range(1, 7):
        expected = (
            n
            * c ** (n - 1)
            * c_prime
            / (math.factorial(n) * math.factorial(n + 1))
        )
        assert jets.parameter_jet(n, 1, eta) == pytest.approx(
            expected, rel=2e-12, abs=1e-300
        )


def test_parameter_jets_satisfy_independent_ftc_compatibility(jets) -> None:
    left = -0.43
    right = 0.58
    # AxisCoefficientSpace.Compatible is exactly the FTC identity between
    # successive stored parameter jets.  Numerical quadrature here is an
    # independent check of the analytic Taylor recurrence used in production.
    for n, m in ((1, 0), (2, 1), (4, 2)):
        integral, _ = quad(
            lambda x: jets.parameter_jet(n, m + 1, x),
            left,
            right,
            epsabs=2e-11,
            epsrel=2e-11,
            limit=100,
        )
        difference = jets.parameter_jet(n, m, right) - jets.parameter_jet(n, m, left)
        assert integral == pytest.approx(difference, rel=2e-9, abs=2e-11)


def test_normalized_coordinate_uses_pinned_axis_weight(jets) -> None:
    n, m, eta = 3, 2, 0.17
    expected_weight = (
        (1.0 / 20.0) ** n
        * jets.epsilon ** (-m)
        * math.factorial(m)
        * math.comb(n + m, m)
        / ((n + 1.0) ** 2 * (m + 1.0) ** 2)
    )
    assert axis_weight(jets.epsilon, n, m) == pytest.approx(expected_weight, rel=2e-15)
    assert jets.normalized_axis_coordinate(n, m, eta) == pytest.approx(
        jets.parameter_jet(n, m, eta) / expected_weight, rel=2e-15
    )


def test_angular_jet_guards_fail_closed(jets) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        jets.parameter_jet(0, 0, 1.100001)
