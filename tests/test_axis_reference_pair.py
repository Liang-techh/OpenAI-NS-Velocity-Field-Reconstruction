import math

import pytest

from openai_ns_reconstruction.axis_reference_pair import (
    actual_schedule_reference_pair,
    radial_divisor,
)
from openai_ns_reconstruction.natural_axis import A, H, L, axis_U, chi, d
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_axis_pressure import (
    axis_pressure,
    axis_pressure_derivative,
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_schedule_reference_uses_certified_sigma_and_stays_fail_closed() -> None:
    reference = actual_schedule_reference_pair(_schedule_data(), 0.05)

    assert reference.sigma == reference.margin_witness.cutoff.sigma
    assert reference.sigma > 0.0
    assert reference.paper_exact is False


def test_reference_angular_coefficients_match_independent_closed_form() -> None:
    reference = actual_schedule_reference_pair(_schedule_data(), 0.05)
    eta = 0.37
    c = chi(reference.data.h, reference.j, reference.sigma, eta)

    # AxisReference.reference_coefficient gives this closed form directly;
    # production code evaluates the recurrence instead, avoiding a test that
    # simply invokes the same implementation twice.
    for n in range(9):
        expected = (-c / 2.0) ** n / (math.factorial(n) * math.factorial(n + 1))
        assert reference.phi_coefficient(n, eta) == pytest.approx(expected, rel=2e-15, abs=0.0)


def test_reference_angular_coefficients_satisfy_resolvent_recurrence() -> None:
    reference = actual_schedule_reference_pair(_schedule_data(), 0.05)
    eta = -0.4
    c = reference.chi0(eta)

    assert reference.phi_coefficient(0, eta) == 1.0
    for n in range(7):
        lhs = radial_divisor(2, n) * reference.phi_coefficient(n + 1, eta)
        rhs = (-c / 2.0) * reference.phi_coefficient(n, eta)
        assert lhs == pytest.approx(rhs, rel=2e-15, abs=1e-300)


def test_reference_axial_coefficient_solves_the_independent_regular_equation() -> None:
    data = _schedule_data()
    j = 0.05
    reference = actual_schedule_reference_pair(data, j)
    eta = 0.23

    # Rebuild NaturalAxisData.Z independently from the actual schedule pressure,
    # rather than calling reference.z_star() as the oracle.
    u_star = axis_U(j, eta)
    pressure = axis_pressure(data, eta)
    pressure_eta = axis_pressure_derivative(data, eta)
    z_star = (
        -A(data.h) * (1.0 - 2.0 * eta * u_star) * u_star
        - 4.0 * H(data.h, j, eta)
        - d(eta) * pressure_eta
        + 4.0 * A(data.h) * eta * pressure
    )
    source = z_star / L(data.h, eta)

    # 2 * (Y u'' + u') = -inverseL*zStar.  For u=u1*Y the left side is
    # 2*radialDivisor(1,0)*u1; all other reference radial coefficients vanish.
    u1 = reference.u_coefficient(1, eta)
    assert 2.0 * radial_divisor(1, 0) * u1 == pytest.approx(-source, rel=2e-13)
    assert reference.u_coefficient(0, eta) == 0.0
    assert reference.u_coefficient(2, eta) == 0.0
    assert reference.u_coefficient(8, eta) == 0.0


def test_reference_axial_data_covers_the_full_pinned_coefficient_window() -> None:
    reference = actual_schedule_reference_pair(_schedule_data(), 0.05)

    # NaturalAxisCoefficients.window is [-11/10,11/10], not merely the original
    # theorem interval [-1,1].  The actual pressure formula is valid here too.
    for eta in (-1.1, 1.1):
        assert math.isfinite(reference.u_coefficient(1, eta))
        assert math.isfinite(reference.phi_coefficient(3, eta))

    with pytest.raises(ValueError, match="pinned window"):
        reference.phi_coefficient(0, 1.100001)


def test_radial_index_and_divisor_guards_fail_closed() -> None:
    reference = actual_schedule_reference_pair(_schedule_data(), 0.05)

    with pytest.raises(ValueError, match="nonnegative integer"):
        reference.phi_coefficient(-1, 0.0)
    with pytest.raises(ValueError, match="positive integer"):
        radial_divisor(0, 0)
