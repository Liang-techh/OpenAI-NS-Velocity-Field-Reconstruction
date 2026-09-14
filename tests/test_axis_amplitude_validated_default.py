from dataclasses import replace
from decimal import (
    Decimal,
    DivisionByZero,
    Inexact,
    InvalidOperation,
    Overflow,
    ROUND_CEILING,
    Underflow,
    localcontext,
)
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_coefficient_amplitude import (
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def crossing():
    amplitude = actual_schedule_amplitude_log_state(_schedule_data(), 0.05)
    eta = -0.02
    phase = amplitude.phase_enclosure(eta)
    enclosure = amplitude.default_log_amplitude_enclosure(eta)
    legacy_log = amplitude.legacy_log_amplitude_diagnostic(eta)
    return amplitude, eta, phase, enclosure, legacy_log


def test_default_phase_crossing_is_enclosed_and_beats_legacy_diagnostic(crossing) -> None:
    amplitude, _, phase, _, legacy_log = crossing
    assert amplitude.phase_absolute_tolerance == Fraction(1, 10**12)
    assert phase.error_bound <= amplitude.phase_absolute_tolerance
    assert float(phase.integral_estimate) == pytest.approx(
        0.05073499503328993,
        rel=0,
        abs=1e-10,
    )

    # Recover the old 4001-point phase from its affine log-amplitude output.
    with localcontext() as context:
        context.prec = 96
        legacy_phase = (legacy_log + amplitude.C_exponent) / amplitude.Lambda
    assert float(legacy_phase) == pytest.approx(
        1.521550417558555,
        rel=0,
        abs=1e-12,
    )
    assert not phase.lower <= Fraction(legacy_phase) <= phase.upper


def test_default_log_midpoint_and_scaled_interval_are_explicit(crossing) -> None:
    amplitude, eta, phase, enclosure, _ = crossing
    point = amplitude.log_amplitude(eta)
    lower, upper = enclosure.decimal_bounds()

    assert enclosure.paper_exact is False
    assert amplitude.phase_has_error_enclosure is True
    assert enclosure.phase_lower == phase.lower
    assert enclosure.phase_upper == phase.upper
    assert enclosure.width == 2 * Fraction(amplitude.Lambda) * phase.error_bound
    assert enclosure.width > Fraction(1, 10**12)
    assert Fraction(lower) <= Fraction(point) <= Fraction(upper)
    assert point == enclosure.midpoint_decimal()

    origin = amplitude.default_log_amplitude_enclosure(0.0)
    assert origin.lower == origin.upper == -Fraction(amplitude.C_exponent)
    assert amplitude.log_amplitude(0.0) == amplitude.C_exponent.copy_negate()


def test_default_phase_tolerance_and_hostile_decimal_context(crossing) -> None:
    amplitude, eta, phase, enclosure, _ = crossing
    with pytest.raises(TypeError, match="phase_absolute_tolerance"):
        replace(amplitude, phase_absolute_tolerance=1e-12)
    with pytest.raises(ValueError, match="phase_absolute_tolerance"):
        replace(amplitude, phase_absolute_tolerance=Fraction(0))

    baseline_point = amplitude.log_amplitude(eta)
    baseline_interval = amplitude.default_log_amplitude_enclosure(eta)
    with localcontext() as context:
        context.prec = 3
        context.rounding = ROUND_CEILING
        context.Emax = 2
        context.Emin = -2
        for signal in (
            DivisionByZero,
            InvalidOperation,
            Overflow,
            Underflow,
            Inexact,
        ):
            context.traps[signal] = True
        # The phase was evaluated above; this call checks the cached exact
        # result and the fresh Decimal contexts used by the affine midpoint.
        assert amplitude.phase_enclosure(eta) == phase
        assert amplitude.default_log_amplitude_enclosure(eta) == baseline_interval
        assert amplitude.log_amplitude(eta) == baseline_point

