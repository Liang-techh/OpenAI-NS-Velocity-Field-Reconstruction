from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_amplitude_log_scale import AmplitudeLogSource
from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_mixed_scale import (
    MixedScaleCoefficient,
)
from openai_ns_reconstruction.axis_coefficient_profile_prefix import (
    FormalAxisProfilePrefix,
    formal_axis_profile_prefix,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_phase_log_enclosure import (
    RationalLogAmplitudeEnclosure,
)
from openai_ns_reconstruction.natural_axis_wide import wide_natural_profile_prefix
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _source() -> AmplitudeLogSource:
    return AmplitudeLogSource(
        h=Fraction(1, 100),
        j=Fraction(1, 20),
        sigma=Fraction(1, 10),
        eta=Fraction(0),
        enclosure=RationalLogAmplitudeEnclosure(
            phase_lower=Fraction(0),
            phase_upper=Fraction(0),
            Lambda=Decimal("3"),
            log_C=Decimal(0),
        ),
        midpoint=Decimal("0.5"),
    )


def _formal(source: AmplitudeLogSource, *, eta: float = 0.0) -> FormalAxisProfilePrefix:
    coefficient = MixedScaleCoefficient(
        {(q, 0): Decimal(q + 1) for q in (0, 1, 2, 4, 9)}
    )
    return FormalAxisProfilePrefix(
        angular=coefficient,
        axial=coefficient,
        axial_average=coefficient,
        pressure=coefficient,
        Lambda=source.enclosure.Lambda,
        amplitude_log=source.midpoint,
        max_n=0,
        radial_order=0,
        eta_order=0,
        Y=Decimal(0),
        eta=eta,
        amplitude_log_source=source,
    )


def test_synthetic_formal_terms_attach_one_source_for_all_q_channels() -> None:
    source = _source()
    terms = _formal(source).angular_terms_log()

    assert set(terms) == {(0, 0), (1, 0), (2, 0), (4, 0), (9, 0)}
    for (q, _), term in terms.items():
        metadata = term.amplitude_log_scale
        assert metadata is not None
        assert metadata.source is source
        assert metadata.q == q
        assert metadata.returned_log_scale == term.log_scale
        expected = Decimal(0) if q == 0 else Decimal(q) * source.midpoint
        assert term.log_scale == expected


def test_formal_profile_rejects_source_with_different_exact_eta() -> None:
    with pytest.raises(ValueError, match="eta"):
        _formal(_source(), eta=0.25)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def solver():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return formal_axis_coefficient_solver(x1)


def test_actual_formal_and_wide_factories_propagate_source_at_axis(solver) -> None:
    formal = formal_axis_profile_prefix(solver, 0, Decimal(0), 0.0)
    source = formal.amplitude_log_source
    assert source is not None
    assert source.midpoint == formal.amplitude_log
    assert source.enclosure.Lambda == formal.Lambda

    formal_terms = formal.angular_terms_log()
    assert formal_terms
    assert all(term.amplitude_log_scale is not None for term in formal_terms.values())

    wide = wide_natural_profile_prefix(solver, 0, Decimal(0), 0.0)
    wide_source = wide.amplitude_log_source
    assert wide_source is not None
    source.assert_compatible(wide_source)
    wide_terms = wide.terms_log("U")
    assert wide_terms
    assert all(term.amplitude_log_scale is not None for term in wide_terms.values())
