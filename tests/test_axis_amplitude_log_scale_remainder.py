from dataclasses import replace
from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_amplitude_log_scale import (
    AmplitudeLogSource,
    AmplitudePowerLogScale,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_remainder import (
    MixedScaleFirstPicardNaturalRemainderCoefficientJet,
    SecondPicardCoefficientJet,
    wide_first_picard_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.axis_phase_log_enclosure import RationalLogAmplitudeEnclosure


LAMBDA = Decimal("10")
AMPLITUDE_LOG = Decimal("-1.5")


def _source() -> AmplitudeLogSource:
    return AmplitudeLogSource(
        h=Fraction(1, 100),
        j=Fraction(1, 10),
        sigma=Fraction(1),
        eta=Fraction(0),
        enclosure=RationalLogAmplitudeEnclosure(
            phase_lower=Fraction(-3, 20),
            phase_upper=Fraction(-149, 1000),
            Lambda=LAMBDA,
            log_C=Decimal("0"),
        ),
        midpoint=AMPLITUDE_LOG,
    )


def _remainder_jet(source=None):
    return MixedScaleFirstPicardNaturalRemainderCoefficientJet(
        ordinary=(Decimal("1"),) * 6,
        pressure_linear=(Decimal("2"), Decimal("-3"), Decimal("4"), Decimal("5"), Decimal("6")),
        pressure_square_inverse_lambda_cubed=Decimal("-7"),
        Lambda=LAMBDA,
        amplitude_log=AMPLITUDE_LOG,
        amplitude_log_source=source,
    )


def _second_jet(source=None):
    return SecondPicardCoefficientJet(
        ordinary=(Decimal("1"),) * 7,
        pressure_linear=(Decimal("0"), Decimal("2"), Decimal("-3"), Decimal("4"), Decimal("5"), Decimal("6")),
        pressure_square_inverse_lambda_fourth=Decimal("-7"),
        Lambda=LAMBDA,
        amplitude_log=AMPLITUDE_LOG,
        amplitude_log_source=source,
    )


def _assert_power_metadata(term, source, q):
    assert isinstance(term.amplitude_log_scale, AmplitudePowerLogScale)
    assert term.amplitude_log_scale.source is source
    assert term.amplitude_log_scale.q == q
    assert term.amplitude_log_scale.returned_log_scale == term.log_scale


def test_synthetic_remainder_and_second_picard_attach_q2_q4_or_legacy_none():
    source = _source()
    remainder = _remainder_jet(source)
    second = _second_jet(source)

    for term in remainder.pressure_linear_terms_log():
        _assert_power_metadata(term, source, 2)
    _assert_power_metadata(remainder.pressure_square_term_log(), source, 4)

    for term in second.pressure_linear_terms_log():
        if term.sign:
            _assert_power_metadata(term, source, 2)
        else:
            assert term.amplitude_log_scale is None
    _assert_power_metadata(second.pressure_square_term_log(), source, 4)

    legacy_remainder = _remainder_jet()
    legacy_second = _second_jet()
    assert all(term.amplitude_log_scale is None for term in legacy_remainder.pressure_linear_terms_log())
    assert legacy_remainder.pressure_square_term_log().amplitude_log_scale is None
    assert all(term.amplitude_log_scale is None for term in legacy_second.pressure_linear_terms_log())
    assert legacy_second.pressure_square_term_log().amplitude_log_scale is None


@pytest.fixture(scope="module")
def actual_remainder():
    data = TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )
    x1 = actual_schedule_wide_first_picard_state(data, 0.05)
    return wide_first_picard_remainder_state(x1)


def test_actual_eta_zero_jet_pairs_share_the_amplitude_source(actual_remainder):
    source = actual_remainder.x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    first_angular, first_axial = actual_remainder.jet_pair(0, 0, 0.0)
    second_angular, second_axial = actual_remainder.second_picard_jet_pair(0, 0, 0.0)

    for jet in (first_angular, first_axial, second_angular, second_axial):
        jet.amplitude_log_source.assert_compatible(source)
        assert jet.amplitude_log == source.midpoint


def test_validate_metadata_rejects_a_branch_from_a_different_eta(actual_remainder):
    source = actual_remainder.amplitude_log_source(0.0)
    actual_log = actual_remainder.amplitude_log(0.0)
    wrong_eta_source = replace(source, eta=Fraction(1, 10))
    branch = replace(
        _remainder_jet(),
        Lambda=actual_remainder.Lambda,
        amplitude_log=actual_log,
        amplitude_log_source=wrong_eta_source,
    )

    with pytest.raises(ValueError, match="eta mismatch"):
        actual_remainder._validate_metadata(
            branch,
            name="synthetic branch",
            Lambda=actual_remainder.Lambda,
            amplitude_log=actual_log,
            amplitude_log_source=source,
        )
