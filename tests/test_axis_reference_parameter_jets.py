import math

import pytest

from openai_ns_reconstruction.axis_reference_parameter_jets import (
    actual_schedule_angular_reference_jets,
    coefficient_weight,
)
from openai_ns_reconstruction.natural_axis import H
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_zeroth_parameter_jet_is_the_landed_reference_coefficient() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)

    for eta in (-1.1, -0.37, 0.0, 0.41, 1.1):
        for n in range(8):
            assert jets.parameter_jet(n, 0, eta) == pytest.approx(
                jets.reference.phi_coefficient(n, eta),
                rel=3e-14,
                abs=1e-300,
            )


def test_first_parameter_jet_matches_independent_closed_formula() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)
    reference = jets.reference
    eta = 0.37

    h_value = H(reference.data.h, reference.j, eta)
    # H = j + (D+4)eta - j eta^2 - 4 eta^3.
    h_prime = (0.5 - reference.data.h + 4.0) - 2.0 * reference.j * eta - 12.0 * eta * eta
    sigma2 = reference.sigma * reference.sigma
    chi_prime = 2.0 * h_value * h_prime * sigma2 / (h_value * h_value + sigma2) ** 2

    # n=1 gives phi0[1] = -chi/4 exactly.
    assert jets.parameter_jet(1, 1, eta) == pytest.approx(-chi_prime / 4.0, rel=3e-13)


def test_constant_reference_mode_has_exact_zero_higher_parameter_jets() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)

    assert jets.parameter_jet(0, 0, -0.8) == 1.0
    for order in range(1, 9):
        assert jets.parameter_jet(0, order, -0.8) == 0.0


def test_high_parameter_jets_are_analytic_not_default_zero() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)

    values = [jets.parameter_jet(3, order, 0.23) for order in range(1, 7)]
    assert all(math.isfinite(value) for value in values)
    assert any(value != 0.0 for value in values[2:])


def test_normalized_coordinates_use_pinned_axis_weight() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)
    n, m, eta = 2, 3, -0.31

    weight = coefficient_weight(jets.epsilon, n, m)
    assert weight > 0.0
    assert jets.normalized_jet(n, m, eta) * weight == pytest.approx(
        jets.parameter_jet(n, m, eta),
        rel=3e-15,
    )
    assert coefficient_weight(jets.epsilon, 0, 0) == 1.0


def test_factory_uses_actual_schedule_radius_and_stays_fail_closed() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)

    assert jets.epsilon > 0.0
    assert jets.paper_exact is False
    assert jets.coefficient_space_membership_certified is False


def test_invalid_indices_and_window_fail_closed() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)

    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(1, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        jets.parameter_jet(1, 1, 1.100001)
