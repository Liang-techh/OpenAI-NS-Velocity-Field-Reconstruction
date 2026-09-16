from dataclasses import replace
from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_amplitude_log_scale import (
    AmplitudeLogSource,
    AmplitudePowerLogScale,
)
from openai_ns_reconstruction.axis_coefficient_amplitude import (
    SignedLogCoefficientJet,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_angular_nonlinear import (
    MixedScaleFirstPicardAngularNonlinearCoefficientJet,
    wide_first_picard_angular_nonlinear_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_axial_square import (
    MixedScaleFirstPicardAxialSquareCoefficientJet,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_lin2 import (
    MixedScaleFirstPicardLin2CoefficientJet,
)
from openai_ns_reconstruction.axis_phase_log_enclosure import (
    RationalLogAmplitudeEnclosure,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


LAMBDA = Decimal("10")
AMPLITUDE_LOG = Decimal("-1.5")


@pytest.fixture(scope="module")
def amplitude_source() -> AmplitudeLogSource:
    enclosure = RationalLogAmplitudeEnclosure(
        Fraction(0),
        Fraction(1, 100),
        LAMBDA,
        Decimal("2"),
    )
    return AmplitudeLogSource(
        h=Fraction(1, 100),
        j=Fraction(1, 10),
        sigma=Fraction(1),
        eta=Fraction(0),
        enclosure=enclosure,
        midpoint=AMPLITUDE_LOG,
    )


def _assert_power(term, source: AmplitudeLogSource, q: int) -> None:
    assert term.sign != 0
    metadata = term.amplitude_log_scale
    assert isinstance(metadata, AmplitudePowerLogScale)
    assert metadata.q == q
    assert metadata.source is source
    assert metadata.returned_log_scale == term.log_scale


def _common_fields(source: AmplitudeLogSource) -> dict[str, object]:
    return {
        "ordinary_reference": Decimal("1"),
        "ordinary_inverse_lambda_numerator": Decimal("2"),
        "ordinary_inverse_lambda_squared_numerator": Decimal("3"),
        "ordinary_inverse_lambda_cubed_numerator": Decimal("4"),
        "ordinary_inverse_lambda_fourth_numerator": Decimal("5"),
        "pressure_linear_inverse_lambda_numerator": Decimal("6"),
        "pressure_linear_inverse_lambda_squared_numerator": Decimal("-7"),
        "pressure_linear_inverse_lambda_cubed_numerator": Decimal("8"),
        "pressure_square_inverse_lambda_squared_numerator": Decimal("9"),
        "Lambda": LAMBDA,
        "amplitude_log": AMPLITUDE_LOG,
        "amplitude_log_source": source,
    }


def test_axial_square_synthetic_q2_q4_metadata(amplitude_source) -> None:
    jet = MixedScaleFirstPicardAxialSquareCoefficientJet(**_common_fields(amplitude_source))
    for term in jet.pressure_linear_terms_log():
        _assert_power(term, amplitude_source, 2)
    _assert_power(jet.pressure_square_term_log(), amplitude_source, 4)


def test_lin2_synthetic_q2_metadata_and_legacy_none(amplitude_source) -> None:
    fields = _common_fields(amplitude_source)
    for name in (
        "ordinary_inverse_lambda_cubed_numerator",
        "ordinary_inverse_lambda_fourth_numerator",
        "pressure_linear_inverse_lambda_squared_numerator",
        "pressure_linear_inverse_lambda_cubed_numerator",
        "pressure_square_inverse_lambda_squared_numerator",
    ):
        fields.pop(name)
    jet = MixedScaleFirstPicardLin2CoefficientJet(**fields)
    _assert_power(jet.pressure_linear_term_log(), amplitude_source, 2)

    legacy = MixedScaleFirstPicardLin2CoefficientJet(
        **{**fields, "amplitude_log_source": None}
    )
    assert legacy.pressure_linear_term_log().amplitude_log_scale is None


def test_angular_nonlinear_synthetic_q2_metadata(amplitude_source) -> None:
    jet = MixedScaleFirstPicardAngularNonlinearCoefficientJet(
        ordinary=(Decimal(1), Decimal(2), Decimal(3), Decimal(4), Decimal(5)),
        pressure_linear=(Decimal(6), Decimal(-7), Decimal(8)),
        Lambda=LAMBDA,
        amplitude_log=AMPLITUDE_LOG,
        amplitude_log_source=amplitude_source,
    )
    for term in jet.pressure_linear_terms_log():
        _assert_power(term, amplitude_source, 2)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)


def test_actual_eta_zero_factories_attach_one_source(x1) -> None:
    expected = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    from openai_ns_reconstruction.axis_coefficient_wide_first_picard_angular_nonlinear import (
        wide_first_picard_angular_nonlinear_state,
    )
    from openai_ns_reconstruction.axis_coefficient_wide_first_picard_axial_square import (
        wide_first_picard_axial_square_state,
    )
    from openai_ns_reconstruction.axis_coefficient_wide_first_picard_lin2 import (
        wide_first_picard_lin2_state,
    )

    jets = (
        wide_first_picard_axial_square_state(x1).jet(0, 0, 0.0),
        wide_first_picard_lin2_state(x1).jet(0, 0, 0.0),
        wide_first_picard_angular_nonlinear_state(x1).quad1_jet(0, 0, 0.0),
    )
    for jet in jets:
        assert jet.amplitude_log_source is not None
        jet.amplitude_log_source.assert_compatible(expected)
        assert jet.amplitude_log == expected.midpoint


def test_actual_angular_pressure_identity_rejects_wrong_source(x1) -> None:
    nonlinear = wide_first_picard_angular_nonlinear_state(x1)
    row = next(
        n
        for n in range(1, 9)
        if x1.jet_pair(n, 0, 0.0)[1].pressure_over_two_lambda.sign != 0
    )
    angular, axial = x1.jet_pair(row, 0, 0.0)
    pressure = axial.pressure_over_two_lambda
    assert pressure.amplitude_log_scale is not None
    expected = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    wrong_source = replace(expected, eta=Fraction(1, 100))
    wrong_metadata = wrong_source.power(2, pressure.log_scale)
    wrong_pressure = SignedLogCoefficientJet(
        sign=pressure.sign,
        log_scale=pressure.log_scale,
        log_factor=pressure.log_factor,
        amplitude_log_scale=wrong_metadata,
    )
    wrong_axial = replace(axial, pressure_over_two_lambda=wrong_pressure)
    with pytest.raises(ValueError, match="amplitude log source"):
        nonlinear._validate_jet_pair_identity(angular, wrong_axial, 0.0)
