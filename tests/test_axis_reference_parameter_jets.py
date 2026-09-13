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
    h_prime = (0.5 - reference.data.h + 4.0) - 2.0 * reference.j * eta - 12.0 * eta * eta
    sigma2 = reference.sigma * reference.sigma
    chi_prime = 2.0 * h_value * h_prime * sigma2 / (h_value * h_value + sigma2) ** 2

    # n=1 gives phi0[1] = -chi/4 exactly.
    assert jets.parameter_jet(1, 1, eta) == pytest.approx(-chi_prime / 4.0, rel=3e-13)


def test_constant_reference_mode_has_exact_zero_higher_parameter_jets() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)

    assert jets.parameter_jet(0, 0, -0.8) == 1.0
    for order in range(1, 9):
        assert jets.parameter_jet_decimal(0, order, -0.8) == 0
        assert jets.parameter_jet(0, order, -0.8) == 0.0


def test_high_parameter_jets_survive_binary64_chi_cancellation() -> None:
    jets = actual_schedule_angular_reference_jets(_schedule_data(), 0.05)
    eta = 0.23

    # On this conservative schedule the direct binary64 chi evaluation rounds
    # to one away from the H-root; the derivative family must not inherit that
    # rounding artefact as fake zero higher jets.
    assert jets.reference.chi0(eta) == 1.0
    wide_values = [jets.parameter_jet_decimal(3, order, eta) for order in range(3, 7)]
    values = [jets.parameter_jet(3, order, eta) for order in range(3, 7)]

    assert all(value.is_finite() and value != 0 for value in wide_values)
    assert all(math.isfinite(value) and value != 0.0 for value in values)
    assert values == [float(value) for value in wide_values]


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
