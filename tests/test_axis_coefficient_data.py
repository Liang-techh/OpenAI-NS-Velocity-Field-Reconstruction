import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_data import (
    LEAN_AXIS_DATA_FIELD_NAMES,
    actual_schedule_axis_coefficient_data,
)
from openai_ns_reconstruction.axis_coefficient_operators import actual_schedule_coefficient_operators
from openai_ns_reconstruction.axis_coefficient_reference_state import actual_schedule_reference_axis_state
from openai_ns_reconstruction.natural_axis import A, D, H, L, W, axis_U, d, real_gradient
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_axis_pressure import axis_pressure, axis_pressure_derivative


ROOT = Path(__file__).resolve().parents[1]


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def reference():
    return actual_schedule_reference_axis_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def coeff_data(reference):
    return actual_schedule_axis_coefficient_data(reference)


def _z_star_value(data: TailData, j: float, eta: float) -> float:
    pressure = axis_pressure(data, eta)
    pressure_eta = axis_pressure_derivative(data, eta)
    u = axis_U(j, eta)
    return (
        -A(data.h) * (1.0 - 2.0 * eta * u) * u
        - 4.0 * H(data.h, j, eta)
        - d(eta) * pressure_eta
        + 4.0 * A(data.h) * eta * pressure
    )


def test_axis_data_layout_and_theorem_selected_scalars(reference, coeff_data) -> None:
    assert coeff_data.lean_field_names == LEAN_AXIS_DATA_FIELD_NAMES
    assert LEAN_AXIS_DATA_FIELD_NAMES == (
        "A",
        "D",
        "h",
        "one",
        "eta",
        "d",
        "inverseL",
        "uStar",
        "uStarEta",
        "wStar",
        "hStar",
        "normalizedGradient",
        "zStar",
    )
    assert coeff_data.A == pytest.approx(A(reference.reference.data.h))
    assert coeff_data.D == pytest.approx(D(reference.reference.data.h))
    assert coeff_data.h == pytest.approx(reference.reference.data.h)
    assert coeff_data.j == pytest.approx(reference.reference.j)
    assert coeff_data.sigma == pytest.approx(reference.reference.sigma)
    assert coeff_data.epsilon == reference.epsilon
    assert coeff_data.paper_exact is False
    assert coeff_data.global_axis_norm_certified is False


def test_fixed_fields_are_radially_constant_and_match_official_real_fields(reference, coeff_data) -> None:
    eta = -0.23
    data = reference.reference.data
    j = reference.reference.j
    expected = {
        "one": 1.0,
        "eta": eta,
        "d": d(eta),
        "inverseL": 1.0 / L(data.h, eta),
        "uStar": axis_U(j, eta),
        "uStarEta": 4.0,
        "wStar": W(data.h, j, eta),
        "hStar": H(data.h, j, eta),
        "normalizedGradient": real_gradient(data.h, j, reference.reference.sigma, eta),
        "zStar": _z_star_value(data, j, eta),
    }
    for name, value in expected.items():
        state = getattr(coeff_data, name)
        assert state.epsilon == reference.epsilon
        assert state.jet(0, 0, eta) == pytest.approx(value, rel=3e-12, abs=3e-11)
        for n in (1, 2, 5):
            for m in (0, 1, 2):
                assert state.jet(n, m, eta) == 0.0


def test_parameter_jets_match_independent_centered_differences(reference, coeff_data) -> None:
    eta = 0.19
    step = 2.0e-5
    data = reference.reference.data
    j = reference.reference.j
    sigma = reference.reference.sigma

    point_functions = {
        "inverseL": lambda x: 1.0 / L(data.h, x),
        "wStar": lambda x: W(data.h, j, x),
        "hStar": lambda x: H(data.h, j, x),
        "normalizedGradient": lambda x: real_gradient(data.h, j, sigma, x),
        "zStar": lambda x: _z_star_value(data, j, x),
    }
    for name, fn in point_functions.items():
        finite_difference = (fn(eta + step) - fn(eta - step)) / (2.0 * step)
        assert getattr(coeff_data, name).jet(0, 1, eta) == pytest.approx(
            finite_difference,
            rel=3e-6,
            abs=2e-6,
        )


def test_landed_coefficient_operators_accept_every_axis_data_field(reference, coeff_data) -> None:
    operators = actual_schedule_coefficient_operators(reference)
    eta = -0.11
    for name in LEAN_AXIS_DATA_FIELD_NAMES[3:]:
        state = getattr(coeff_data, name)
        averaged = operators.average(state)
        assert averaged.epsilon == reference.epsilon
        assert averaged.jet(0, 0, eta) == pytest.approx(state.jet(0, 0, eta))


def test_factory_rejects_non_reference_objects() -> None:
    with pytest.raises(TypeError, match="ActualScheduleReferenceAxisState"):
        actual_schedule_axis_coefficient_data(object())


def test_axis_data_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_data.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisData" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
