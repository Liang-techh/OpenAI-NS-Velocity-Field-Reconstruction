from decimal import Decimal, localcontext
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_angular_square import (
    wide_first_picard_angular_square_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_source import (
    ActualScheduleWideFirstPicardNaturalSourceState,
    MixedScaleFirstPicardNaturalSourceCoefficientJet,
    wide_first_picard_natural_source_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


DECIMAL_PRECISION = 96
FIELDS = (
    "reference_factor",
    "inverse_lambda_numerator",
    "inverse_lambda_squared_numerator",
    "inverse_lambda_cubed_numerator",
    "inverse_lambda_fourth_numerator",
)
SQUARE_FIELDS = (
    "reference",
    "inverse_lambda_numerator",
    "inverse_lambda_squared_numerator",
    "inverse_lambda_cubed_numerator",
    "inverse_lambda_fourth_numerator",
)


def _schedule_data() -> TailData:
    return TailData(OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0), h=0.01)


@pytest.fixture(scope="module")
def source_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_natural_source_state(x1)


def test_source_x1_is_bound_to_the_genuine_first_picard_state(source_state) -> None:
    assert isinstance(source_state, ActualScheduleWideFirstPicardNaturalSourceState)
    assert source_state.x1.picard_x1_materialized is True
    assert source_state.natural_source_x1_materialized is True
    assert source_state.natural_remainder_x1_materialized is False
    assert source_state.picard_x2_materialized is False
    assert source_state.global_axis_norm_certified is False
    assert source_state.paper_exact is False
    assert source_state.Lambda == source_state.amplitude.Lambda
    assert source_state.epsilon == source_state.amplitude.epsilon


def test_zeroth_eta_derivative_is_exactly_a2_times_landed_phi1_square(source_state) -> None:
    n, eta = 4, 0.03
    source = source_state.jet(n, 0, eta)
    square = wide_first_picard_angular_square_state(source_state.x1).jet(n, 0, eta)
    assert isinstance(source, MixedScaleFirstPicardNaturalSourceCoefficientJet)
    for source_name, square_name in zip(FIELDS, SQUARE_FIELDS):
        assert getattr(source, source_name) == getattr(square, square_name)
    assert source.amplitude_log == source_state.amplitude.log_amplitude(eta)


def test_first_eta_derivative_obeys_exact_product_rule_for_a2(source_state) -> None:
    n, eta = 3, -0.04
    source = source_state.jet(n, 1, eta)
    square_state = wide_first_picard_angular_square_state(source_state.x1)
    square0 = square_state.jet(n, 0, eta)
    square1 = square_state.jet(n, 1, eta)
    gradient = source_state.amplitude.data.normalizedGradient.jet(0, 0, eta)
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        bell1 = +(
            Decimal(2)
            * source_state.Lambda
            * Decimal.from_float(float(gradient))
        )
        for source_name, square_name in zip(FIELDS, SQUARE_FIELDS):
            expected = +(
                getattr(square1, square_name)
                + bell1 * getattr(square0, square_name)
            )
            assert getattr(source, source_name) == expected


def test_source_keeps_theorem_amplitude_and_inverse_lambda_scales_split(source_state) -> None:
    selected = next(
        (
            source_state.jet(n, m, 0.07)
            for n in range(1, 10)
            for m in range(0, 3)
            if any(value != 0 for value in source_state.jet(n, m, 0.07).normalized_factors_decimal())
        ),
        None,
    )
    assert selected is not None
    terms = selected.source_terms_log()
    assert len(terms) == 5
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        expected_scale = +(Decimal(2) * selected.amplitude_log)
    for power, (factor, logged) in enumerate(zip(selected.normalized_factors_decimal(), terms)):
        assert logged.sign == (0 if factor == 0 else (1 if factor > 0 else -1))
        if factor == 0:
            assert logged.log_scale is None
            assert logged.log_factor is None
            continue
        assert logged.log_scale == expected_scale
        with localcontext() as ctx:
            ctx.prec = DECIMAL_PRECISION
            expected_factor = +(
                abs(factor).ln()
                - Decimal(power) * selected.Lambda.ln()
            )
        assert logged.log_factor == expected_factor


def test_source_rejects_surrogates_and_invalid_coordinates(source_state) -> None:
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_natural_source_state(object())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        source_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        source_state.jet(1, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        source_state.jet(1, 0, 2.0)


def test_source_provenance_is_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_source.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    truth = layer["truth_boundary"]
    assert truth["first_picard_x1_materialized"] is True
    assert truth["natural_source_x1_materialized"] is True
    assert truth["natural_remainder_x1_materialized"] is False
    assert truth["picard_x2_materialized"] is False
    assert truth["fixed_point_materialized"] is False
    assert truth["global_axis_norm_certified"] is False
    assert "pressure" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
