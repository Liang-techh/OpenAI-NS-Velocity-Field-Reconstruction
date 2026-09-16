from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_amplitude_log_scale import (
    AmplitudePowerLogScale,
    AmplitudeLogSource,
)
from openai_ns_reconstruction.axis_phase_log_enclosure import (
    RationalLogAmplitudeEnclosure,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_dot import (
    MixedScaleFirstPicardSlow2AverageDotCoefficientJet,
    wide_first_picard_slow2_average_dot_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_mixed import (
    MixedScaleFirstPicardSlow2AverageMixedCoefficientJet,
    wide_first_picard_slow2_average_mixed_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_axial_quadratic import (
    MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet,
    wide_first_picard_slow2_axial_quadratic_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_param import (
    MixedScaleFirstPicardSlow2ParamCoefficientJet,
    wide_first_picard_slow2_param_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


LAMBDA = Decimal("10")
AMPLITUDE_LOG = Decimal("-1.5")
BRANCHES = (
    MixedScaleFirstPicardSlow2ParamCoefficientJet,
    MixedScaleFirstPicardSlow2AverageDotCoefficientJet,
    MixedScaleFirstPicardSlow2AverageMixedCoefficientJet,
    MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet,
)


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


def _fields(source: AmplitudeLogSource) -> dict[str, Decimal | AmplitudeLogSource]:
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


def _assert_term_scale(term, source: AmplitudeLogSource, power: int) -> None:
    assert term.sign != 0
    metadata = term.amplitude_log_scale
    assert isinstance(metadata, AmplitudePowerLogScale)
    assert metadata.q == power
    assert metadata.source is source
    assert metadata.returned_log_scale == term.log_scale


@pytest.mark.parametrize("branch", BRANCHES)
def test_synthetic_slow2_branch_terms_keep_q2_q4_source(branch, amplitude_source) -> None:
    jet = branch(**_fields(amplitude_source))
    linear_terms = jet.pressure_linear_terms_log()
    assert len(linear_terms) == 3
    for term in linear_terms:
        _assert_term_scale(term, amplitude_source, 2)
    _assert_term_scale(jet.pressure_square_term_log(), amplitude_source, 4)

    legacy = branch(
        **{**_fields(amplitude_source), "amplitude_log_source": None}
    )
    assert all(term.amplitude_log_scale is None for term in legacy.pressure_linear_terms_log())
    assert legacy.pressure_square_term_log().amplitude_log_scale is None


@pytest.mark.parametrize("branch", BRANCHES)
def test_synthetic_slow2_branch_rejects_source_mismatch(branch, amplitude_source) -> None:
    midpoint_mismatch = _fields(amplitude_source)
    midpoint_mismatch["amplitude_log"] = Decimal("-1")
    with pytest.raises(ValueError, match="midpoint"):
        branch(**midpoint_mismatch)

    lambda_mismatch = _fields(amplitude_source)
    lambda_mismatch["Lambda"] = Decimal("11")
    with pytest.raises(ValueError, match="Lambda"):
        branch(**lambda_mismatch)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)


def test_actual_eta_zero_row_zero_param_and_dot_attach_source(x1) -> None:
    expected = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    for state in (
        wide_first_picard_slow2_param_state(x1),
        wide_first_picard_slow2_average_dot_state(x1),
    ):
        jet = state.jet(0, 0, 0.0)
        assert jet.amplitude_log_source is not None
        jet.amplitude_log_source.assert_compatible(expected)
        assert jet.amplitude_log == expected.midpoint
