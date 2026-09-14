from decimal import Decimal, localcontext
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_source import (
    wide_first_picard_natural_source_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_pressure import (
    ActualScheduleWideFirstPicardNaturalPressureState,
    MixedScaleFirstPicardNaturalPressureCoefficientJet,
    wide_first_picard_natural_pressure_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


DECIMAL_PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def pressure_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    source = wide_first_picard_natural_source_state(x1)
    return wide_first_picard_natural_pressure_state(source)


def test_pressure_x1_is_bound_to_genuine_source_and_axis_data(pressure_state) -> None:
    assert isinstance(
        pressure_state,
        ActualScheduleWideFirstPicardNaturalPressureState,
    )
    assert pressure_state.source.natural_source_x1_materialized is True
    assert pressure_state.natural_pressure_x1_materialized is True
    assert pressure_state.natural_remainder_x1_materialized is False
    assert pressure_state.picard_x2_materialized is False
    assert pressure_state.global_axis_norm_certified is False
    assert pressure_state.paper_exact is False
    assert pressure_state.Lambda == pressure_state.source.Lambda
    assert pressure_state.epsilon == pressure_state.axis_data.epsilon


def test_pinned_radial_maps_are_exact_for_each_inverse_lambda_power(pressure_state) -> None:
    n, m, eta = 4, 1, 0.03
    for power in range(5):
        source = pressure_state.source_factor(power, n - 1, m, eta)
        source_next_eta = pressure_state.source_factor(power, n - 1, m + 1, eta)
        with localcontext() as ctx:
            ctx.prec = DECIMAL_PRECISION
            assert +(pressure_state.primitive_factor(power, n, m, eta) * Decimal(n)) == source
            assert +(
                pressure_state.parameter_primitive_factor(power, n, m, eta)
                * Decimal(n)
            ) == source_next_eta
        assert pressure_state.multiply_y_factor(power, n, m, eta) == source


def test_pressure_input_zeroth_eta_jet_matches_pinned_formula_without_tolerance(
    pressure_state,
) -> None:
    n, eta = 3, -0.04
    d = pressure_state.axis_data
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        A = Decimal.from_float(float(d.A))
        eta0 = Decimal.from_float(float(d.eta.jet(0, 0, eta)))
        d0 = Decimal.from_float(float(d.d.jet(0, 0, eta)))
        for power in range(5):
            primitive = pressure_state.primitive_factor(power, n, 0, eta)
            parameter_primitive = pressure_state.parameter_primitive_factor(
                power, n, 0, eta
            )
            multiply_y = pressure_state.multiply_y_factor(power, n, 0, eta)
            expected = +(
                -Decimal(4) * A * eta0 * primitive
                + d0 * parameter_primitive
                - Decimal(2) * eta0 * multiply_y
            )
            assert pressure_state.pressure_input_factor(power, n, 0, eta) == expected


def test_final_j1_row_identity_is_exact_and_zero_row_is_pinned(pressure_state) -> None:
    m, eta = 2, 0.07
    for power in range(5):
        assert pressure_state.pressure_factor(power, 0, m, eta) == Decimal(0)
        for n in (1, 2, 5):
            with localcontext() as ctx:
                ctx.prec = DECIMAL_PRECISION
                lhs = +(
                    pressure_state.pressure_factor(power, n, m, eta)
                    * Decimal(n * n)
                )
                rhs = pressure_state.pressure_input_factor(
                    power, n - 1, m, eta
                )
            assert lhs == rhs


def test_pressure_keeps_all_nonzero_scales_in_signed_log_form(pressure_state) -> None:
    selected = next(
        (
            pressure_state.jet(n, m, 0.02)
            for n in range(1, 10)
            for m in range(0, 3)
            if any(
                value != 0
                for value in pressure_state.jet(
                    n, m, 0.02
                ).normalized_factors_decimal()
            )
        ),
        None,
    )
    assert selected is not None
    assert isinstance(selected, MixedScaleFirstPicardNaturalPressureCoefficientJet)
    terms = selected.pressure_terms_log()
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        expected_scale = +(Decimal(2) * selected.amplitude_log)
    for power, (factor, logged) in enumerate(
        zip(selected.normalized_factors_decimal(), terms)
    ):
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


def test_pressure_rejects_surrogates_invalid_coordinates_and_scale_indices(
    pressure_state,
) -> None:
    with pytest.raises(
        TypeError,
        match="source must be ActualScheduleWideFirstPicardNaturalSourceState",
    ):
        wide_first_picard_natural_pressure_state(object())
    with pytest.raises(ValueError, match="inverse_lambda_power must lie in 0..4"):
        pressure_state.pressure_factor(5, 1, 0, 0.0)
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        pressure_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        pressure_state.jet(1, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        pressure_state.jet(1, 0, 2.0)


def test_pressure_provenance_is_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_pressure.json"
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
    assert truth["natural_pressure_x1_materialized"] is True
    assert truth["natural_remainder_x1_materialized"] is False
    assert truth["picard_x2_materialized"] is False
    assert truth["fixed_point_materialized"] is False
    assert truth["global_axis_norm_certified"] is False
    assert "remaining angular" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
