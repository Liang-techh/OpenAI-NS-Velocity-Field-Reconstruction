from decimal import Decimal, localcontext

import pytest

from openai_ns_reconstruction.axis_amplitude_log_scale import (
    AmplitudePowerLogScale,
)
from openai_ns_reconstruction.axis_coefficient_amplitude import (
    SignedLogCoefficientJet,
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_axial_remainder import (
    actual_schedule_reference_wide_axial_remainder_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    _pressure_over_two_lambda,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_average import (
    _average_signed_log,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_pressure import (
    _signed_log_term,
)
from openai_ns_reconstruction.axis_coefficient_wide_natural_pressure import (
    actual_schedule_reference_wide_pressure_log_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_natural_source import (
    actual_schedule_reference_natural_source_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def amplitude():
    return actual_schedule_amplitude_log_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def source_state():
    return actual_schedule_reference_natural_source_log_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def pressure_state():
    return actual_schedule_reference_wide_pressure_log_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def axial_state():
    return actual_schedule_reference_wide_axial_remainder_state(_schedule_data(), 0.05)


def _assert_power_metadata(jet, amplitude, q: int) -> None:
    assert jet.sign != 0
    metadata = jet.amplitude_log_scale
    assert isinstance(metadata, AmplitudePowerLogScale)
    assert metadata.q == q
    assert metadata.returned_log_scale == jet.log_scale

    expected_source = amplitude.log_amplitude_source(0.0)
    metadata.source.assert_compatible(expected_source)
    with localcontext() as context:
        context.prec = PRECISION
        expected_scale = +(Decimal(q) * expected_source.midpoint)
    assert jet.log_scale == expected_scale


def _assert_zero(jet) -> None:
    assert jet.sign == 0
    assert jet.log_scale is None
    assert jet.log_factor is None
    assert jet.amplitude_log_scale is None


def test_actual_amplitude_q1_keeps_common_scale_source(amplitude) -> None:
    jet = amplitude.jet_log(0, 0, 0.0)
    _assert_power_metadata(jet, amplitude, 1)


def test_x0_source_q2_keeps_common_scale_source(source_state) -> None:
    jet = source_state.jet_log(0, 0, 0.0)
    _assert_power_metadata(jet, source_state.amplitude, 2)


def test_x0_pressure_q2_has_metadata_and_exact_zero_row(pressure_state) -> None:
    _assert_zero(pressure_state.jet_log(0, 0, 0.0))
    row = next(
        n
        for n in range(1, 9)
        if pressure_state.normalized_factor(n, 0, 0.0) != 0
    )
    _assert_power_metadata(
        pressure_state.jet_log(row, 0, 0.0),
        pressure_state.amplitude,
        2,
    )


def test_x0_axial_pressure_q2_has_metadata_and_exact_zero_row(axial_state) -> None:
    _assert_zero(axial_state.pressure_jet_log(0, 0, 0.0))
    row = next(
        n
        for n in range(1, 9)
        if axial_state.pressure_normalized_factor(n, 0, 0.0) != 0
    )
    _assert_power_metadata(
        axial_state.pressure_jet_log(row, 0, 0.0),
        axial_state.wide_pressure.amplitude,
        2,
    )


def test_metadata_mismatch_and_scale_only_combinators(amplitude) -> None:
    source = amplitude.log_amplitude_source(0.0)
    with localcontext() as context:
        context.prec = PRECISION
        log_scale = +(Decimal(2) * source.midpoint)
    metadata = source.power(2, log_scale)
    bad_metadata = source.power(2, Decimal(0))
    assert bad_metadata.returned_log_scale != log_scale

    with pytest.raises(ValueError, match="returned_log_scale"):
        SignedLogCoefficientJet(
            sign=1,
            log_scale=log_scale,
            log_factor=Decimal(3),
            amplitude_log_scale=bad_metadata,
        )

    jet = SignedLogCoefficientJet(
        sign=1,
        log_scale=log_scale,
        log_factor=Decimal(3),
        amplitude_log_scale=metadata,
    )
    scaled = _pressure_over_two_lambda(jet, amplitude.Lambda)
    averaged = _average_signed_log(jet, 3)
    assert scaled.amplitude_log_scale == metadata
    assert averaged.amplitude_log_scale == metadata
    _assert_zero(_pressure_over_two_lambda(SignedLogCoefficientJet.zero(), amplitude.Lambda))
    _assert_zero(_average_signed_log(SignedLogCoefficientJet.zero(), 3))


def test_x1_pressure_q2_helper_requires_selected_lambda(amplitude) -> None:
    source = amplitude.log_amplitude_source(0.0)
    jet = _signed_log_term(
        Decimal(3),
        amplitude_source=source,
        Lambda=amplitude.Lambda,
        inverse_lambda_power=2,
    )
    _assert_power_metadata(jet, amplitude, 2)

    wrong_lambda = Decimal(1)
    assert wrong_lambda != amplitude.Lambda
    with pytest.raises(ValueError, match="Lambda"):
        _signed_log_term(
            Decimal(3),
            amplitude_source=source,
            Lambda=wrong_lambda,
            inverse_lambda_power=2,
        )
    with pytest.raises(ValueError, match="Lambda"):
        _signed_log_term(
            Decimal(0),
            amplitude_source=source,
            Lambda=wrong_lambda,
            inverse_lambda_power=2,
        )
