from decimal import Decimal
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import actual_schedule_wide_first_picard_state
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2 import (
    ActualScheduleWideFirstPicardSlow2State,
    MixedScaleFirstPicardSlow2CoefficientJet,
    wide_first_picard_slow2_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_axial_quadratic import wide_first_picard_slow2_axial_quadratic_state
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_dot import wide_first_picard_slow2_average_dot_state
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_mixed import wide_first_picard_slow2_average_mixed_state
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_param import wide_first_picard_slow2_param_state
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


FIELDS = (
    "ordinary_reference",
    "ordinary_inverse_lambda_numerator",
    "ordinary_inverse_lambda_squared_numerator",
    "ordinary_inverse_lambda_cubed_numerator",
    "ordinary_inverse_lambda_fourth_numerator",
    "pressure_linear_inverse_lambda_numerator",
    "pressure_linear_inverse_lambda_squared_numerator",
    "pressure_linear_inverse_lambda_cubed_numerator",
    "pressure_square_inverse_lambda_squared_numerator",
)


def _schedule_data() -> TailData:
    return TailData(OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0), h=0.01)


@pytest.fixture(scope="module")
def slow2_state():
    x1 = actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)
    return wide_first_picard_slow2_state(x1)


def _branches(x1):
    return (
        wide_first_picard_slow2_axial_quadratic_state(x1),
        wide_first_picard_slow2_average_dot_state(x1),
        wide_first_picard_slow2_average_mixed_state(x1),
        wide_first_picard_slow2_param_state(x1),
    )


def test_complete_slow2_is_bound_to_one_genuine_x1(slow2_state) -> None:
    assert isinstance(slow2_state, ActualScheduleWideFirstPicardSlow2State)
    assert slow2_state.x1.picard_x1_materialized is True
    assert slow2_state.slow2_all_constituent_branches_materialized is True
    assert slow2_state.slow2_complete is True
    assert slow2_state.natural_remainder_x1_materialized is False
    assert slow2_state.picard_x2_materialized is False
    assert slow2_state.fixed_point_materialized is False
    assert slow2_state.paper_exact is False
    for branch in _branches(slow2_state.x1):
        assert branch.x1 is slow2_state.x1
        assert branch.Lambda == slow2_state.Lambda
        assert branch.epsilon == slow2_state.epsilon


def test_complete_slow2_jet_is_exact_sum_of_four_landed_branches(slow2_state) -> None:
    n, m, eta = 4, 2, 0.03
    actual = slow2_state.jet(n, m, eta)
    constituent_jets = tuple(branch.jet(n, m, eta) for branch in _branches(slow2_state.x1))
    assert isinstance(actual, MixedScaleFirstPicardSlow2CoefficientJet)
    for name in FIELDS:
        assert getattr(actual, name) == sum((getattr(jet, name) for jet in constituent_jets), Decimal(0))
    assert actual.Lambda == slow2_state.Lambda
    assert all(jet.Lambda == actual.Lambda for jet in constituent_jets)
    assert all(jet.amplitude_log == actual.amplitude_log for jet in constituent_jets)


def test_complete_slow2_row_zero_preserves_exact_branch_sum(slow2_state) -> None:
    actual = slow2_state.jet(0, 2, -0.07)
    constituent_jets = tuple(branch.jet(0, 2, -0.07) for branch in _branches(slow2_state.x1))
    for name in FIELDS:
        assert getattr(actual, name) == sum((getattr(jet, name) for jet in constituent_jets), Decimal(0))
        assert getattr(actual, name) == 0


def test_complete_slow2_retains_signed_log_pressure_scales(slow2_state) -> None:
    selected = next(
        (
            slow2_state.jet(n, m, 0.07)
            for n in range(2, 13)
            for m in range(0, 3)
            if any(
                getattr(slow2_state.jet(n, m, 0.07), name) != 0
                for name in (
                    "pressure_linear_inverse_lambda_numerator",
                    "pressure_linear_inverse_lambda_squared_numerator",
                    "pressure_linear_inverse_lambda_cubed_numerator",
                    "pressure_square_inverse_lambda_squared_numerator",
                )
            )
        ),
        None,
    )
    assert selected is not None
    linear = selected.pressure_linear_terms_log()
    for numerator, logged in zip(
        (
            selected.pressure_linear_inverse_lambda_numerator,
            selected.pressure_linear_inverse_lambda_squared_numerator,
            selected.pressure_linear_inverse_lambda_cubed_numerator,
        ),
        linear,
    ):
        assert logged.sign == (0 if numerator == 0 else (1 if numerator > 0 else -1))
        if numerator != 0:
            assert logged.log_scale == Decimal(2) * selected.amplitude_log
    square = selected.pressure_square_term_log()
    numerator = selected.pressure_square_inverse_lambda_squared_numerator
    assert square.sign == (0 if numerator == 0 else (1 if numerator > 0 else -1))
    if numerator != 0:
        assert square.log_scale == Decimal(4) * selected.amplitude_log


def test_complete_slow2_rejects_surrogate_and_delegates_guards(slow2_state) -> None:
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_slow2_state(object())
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        slow2_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        slow2_state.jet(1, 0, 2.0)


def test_complete_slow2_provenance_is_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "references" / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2.json").read_text(encoding="utf-8"))
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    truth = layer["truth_boundary"]
    assert truth["slow2_all_constituent_branches_materialized"] is True
    assert truth["slow2_complete"] is True
    assert truth["natural_remainder_x1_materialized"] is False
    assert truth["picard_x2_materialized"] is False
    assert truth["fixed_point_materialized"] is False
    assert "naturalRemainder(x1)" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
