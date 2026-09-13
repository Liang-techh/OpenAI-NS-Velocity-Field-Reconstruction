from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_coefficient_amplitude import (
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.axis_phase_integral import PhaseIntegrationLimit
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    """Use the same legitimate actual-schedule fixture as amplitude tests."""

    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def amplitude():
    return actual_schedule_amplitude_log_state(_schedule_data(), 0.05)


def test_origin_enclosure_keeps_exact_symbolic_offset(amplitude) -> None:
    tolerance = Fraction(1, 10**20)
    enclosure = amplitude.log_amplitude_enclosure(
        Fraction(0),
        absolute_log_tolerance=tolerance,
    )

    assert enclosure.lower == enclosure.upper == -Fraction(amplitude.C_exponent)
    assert enclosure.width == 0
    assert enclosure.width <= 2 * tolerance
    assert enclosure.paper_exact is False


def test_nonzero_actual_phase_interval_respects_exact_log_tolerance(amplitude) -> None:
    # The actual schedule's Lambda is intentionally very large.  A deliberately
    # wide log tolerance keeps this nonzero phase call within a cheap order-16
    # validated integration while still testing exact Lambda/tolerance coupling.
    tolerance = Fraction(10) ** 777
    enclosure = amplitude.log_amplitude_enclosure(
        Fraction(1, 100),
        absolute_log_tolerance=tolerance,
    )

    assert enclosure.lower < enclosure.upper
    assert enclosure.width <= 2 * tolerance
    lower, upper = enclosure.decimal_bounds()
    assert Fraction(lower) <= enclosure.lower
    assert Fraction(upper) >= enclosure.upper


def test_exact_fraction_inputs_and_integrator_caps_fail_closed(amplitude) -> None:
    tolerance = Fraction(1, 10**100)
    with pytest.raises(TypeError, match="eta must be a Fraction"):
        amplitude.log_amplitude_enclosure(0.0, absolute_log_tolerance=tolerance)
    with pytest.raises(TypeError, match="absolute_log_tolerance must be a Fraction"):
        amplitude.log_amplitude_enclosure(Fraction(0), absolute_log_tolerance=1)
    with pytest.raises(ValueError, match="absolute_log_tolerance must be positive"):
        amplitude.log_amplitude_enclosure(
            Fraction(0),
            absolute_log_tolerance=Fraction(0),
        )

    # The exact phase budget is far below the order-1 cell remainder here;
    # PhaseIntegrationLimit must pass through the amplitude bridge unchanged.
    with pytest.raises(PhaseIntegrationLimit):
        amplitude.log_amplitude_enclosure(
            Fraction(1, 100),
            absolute_log_tolerance=tolerance,
            initial_order=0,
            max_order=1,
            max_cells=8,
            max_depth=8,
        )
