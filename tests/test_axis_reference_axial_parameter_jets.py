import math

import pytest

from openai_ns_reconstruction.axis_reference_axial_parameter_jets import (
    actual_schedule_axial_reference_jets,
)
from openai_ns_reconstruction.axis_reference_parameter_jets import coefficient_weight
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_zeroth_axial_parameter_jet_is_landed_reference_coefficient() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)

    for eta in (-1.1, -0.37, 0.0, 0.41, 1.1):
        assert jets.parameter_jet(1, 0, eta) == jets.reference.u_coefficient(1, eta)


def test_non_degree_one_axial_modes_are_exact_structural_zeros() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)

    for n in (0, 2, 3, 9):
        for m in (0, 1, 4, 8):
            assert jets.parameter_jet(n, m, 0.23) == 0.0
            assert jets.normalized_parameter_taylor_coefficient(n, m, 0.23) == 0.0


def test_first_axial_parameter_jet_matches_independent_centered_oracle() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)
    eta = 0.23
    step = 2.0e-5

    oracle = (
        jets.reference.u_coefficient(1, eta + step)
        - jets.reference.u_coefficient(1, eta - step)
    ) / (2.0 * step)
    value = jets.parameter_jet(1, 1, eta)

    # The production derivative is analytic Taylor algebra.  The centered
    # difference is used only here as an independent regression oracle.
    assert math.isfinite(value)
    assert value == pytest.approx(oracle, rel=3e-6, abs=2e-8)


def test_higher_axial_parameter_jets_are_materialized_from_actual_pressure() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)
    eta = 0.23

    values = [jets.parameter_jet(1, order, eta) for order in (2, 3, 4, 5)]
    assert all(math.isfinite(value) for value in values)
    assert any(value != 0.0 for value in values)


def test_normalized_coordinates_use_same_pinned_axis_weight() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)
    n, m, eta = 1, 3, -0.31

    weight = coefficient_weight(jets.epsilon, n, m)
    assert weight > 0.0
    assert jets.normalized_jet(n, m, eta) * weight == pytest.approx(
        jets.parameter_jet(n, m, eta),
        rel=3e-15,
    )


def test_factory_uses_actual_schedule_radius_and_remains_fail_closed() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)

    assert jets.epsilon > 0.0
    assert jets.paper_exact is False
    assert jets.coefficient_space_membership_certified is False


def test_invalid_indices_and_window_fail_closed() -> None:
    jets = actual_schedule_axial_reference_jets(_schedule_data(), 0.05)

    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(1, -1, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        jets.parameter_jet(True, 0, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        jets.parameter_jet(1, 1, 1.100001)
