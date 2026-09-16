from dataclasses import replace
from decimal import Decimal
from fractions import Fraction
from types import SimpleNamespace

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2 import (
    wide_first_picard_slow2_state,
)
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


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(
        TailData(
            OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
            h=0.01,
        ),
        0.05,
    )


def _zero_jet(x1, source):
    return SimpleNamespace(
        **{name: Decimal(0) for name in FIELDS},
        Lambda=x1.Lambda,
        amplitude_log=source.midpoint,
        amplitude_log_source=source,
    )


def _install_zero_branches(state, x1, sources):
    for attribute, source in sources.items():
        jet = _zero_jet(x1, source)
        branch = SimpleNamespace(
            x1=x1,
            epsilon=x1.epsilon,
            Lambda=x1.Lambda,
            slow2_axial_quadratic_branch_materialized=True,
            slow2_average_dot_branch_materialized=True,
            slow2_average_mixed_branch_materialized=True,
            slow2_param_branch_materialized=True,
            jet=lambda n, m, eta, jet=jet: jet,
        )
        object.__setattr__(state, attribute, branch)


def test_all_four_zero_branches_keep_one_actual_source(x1) -> None:
    state = wide_first_picard_slow2_state(x1)
    source = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    _install_zero_branches(
        state,
        x1,
        {name: source for name in ("axial_quadratic", "average_dot", "average_mixed", "param")},
    )

    result = state.jet(0, 0, 0.0)

    assert result.amplitude_log_source is source
    assert result.amplitude_log == source.midpoint
    assert all(getattr(result, name) == 0 for name in FIELDS)


def test_one_changed_eta_source_is_rejected(x1) -> None:
    state = wide_first_picard_slow2_state(x1)
    source = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    changed = replace(source, eta=Fraction(1, 10))
    _install_zero_branches(
        state,
        x1,
        {
            "axial_quadratic": source,
            "average_dot": changed,
            "average_mixed": source,
            "param": source,
        },
    )

    with pytest.raises(ValueError, match="eta mismatch"):
        state.jet(0, 0, 0.0)


def test_all_consistently_wrong_eta_sources_are_rejected_against_actual_anchor(x1) -> None:
    state = wide_first_picard_slow2_state(x1)
    actual = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    wrong = replace(actual, eta=Fraction(1, 10))
    _install_zero_branches(
        state,
        x1,
        {name: wrong for name in ("axial_quadratic", "average_dot", "average_mixed", "param")},
    )

    with pytest.raises(ValueError, match="eta mismatch"):
        state.jet(0, 0, 0.0)


def test_missing_source_metadata_is_rejected(x1) -> None:
    state = wide_first_picard_slow2_state(x1)
    source = x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(0.0)
    _install_zero_branches(
        state,
        x1,
        {name: source for name in ("axial_quadratic", "average_dot", "average_mixed", "param")},
    )
    missing = state.param.jet(0, 0, 0.0)
    object.__setattr__(missing, "amplitude_log_source", None)

    with pytest.raises(ValueError, match="missing amplitude log source metadata"):
        state.jet(0, 0, 0.0)
