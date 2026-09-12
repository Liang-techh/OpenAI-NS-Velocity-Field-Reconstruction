import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_average import axis_coefficient_average
from openai_ns_reconstruction.axis_coefficient_inverse_dot_product import axis_coefficient_inverse_dot_product
from openai_ns_reconstruction.axis_coefficient_inverse_mixed import axis_coefficient_inverse_mixed
from openai_ns_reconstruction.axis_coefficient_inverse_param_product import axis_coefficient_inverse_param_product
from openai_ns_reconstruction.axis_coefficient_multiply_y import axis_coefficient_multiply_y
from openai_ns_reconstruction.axis_coefficient_operators import (
    LEAN_FIELD_NAMES,
    actual_schedule_coefficient_operators,
)
from openai_ns_reconstruction.axis_coefficient_parameter_primitive import axis_coefficient_parameter_primitive
from openai_ns_reconstruction.axis_coefficient_primitive import axis_coefficient_primitive
from openai_ns_reconstruction.axis_coefficient_product import axis_coefficient_product
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.axis_coefficient_regular_inverse import axis_coefficient_regular_inverse
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


ROOT = Path(__file__).resolve().parents[1]


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def state():
    return actual_schedule_reference_axis_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def operators(state):
    return actual_schedule_coefficient_operators(state)


def _same_jets(actual, expected, *, eta: float = -0.17) -> None:
    for n in range(6):
        for m in range(3):
            assert actual.jet(n, m, eta) == pytest.approx(
                expected.jet(n, m, eta),
                rel=2e-14,
                abs=2e-10,
            )


def test_record_exposes_exact_pinned_lean_field_layout(operators) -> None:
    assert operators.lean_field_names == LEAN_FIELD_NAMES
    assert LEAN_FIELD_NAMES == (
        "product",
        "average",
        "primitive",
        "parameterPrimitive",
        "mulY",
        "j1",
        "j2",
        "param1",
        "param2",
        "dot1",
        "dot2",
        "mixed1",
        "mixed2",
    )
    for name in LEAN_FIELD_NAMES:
        assert callable(getattr(operators, name))


def test_unary_fields_are_exact_assemblies_of_landed_primitives(state, operators) -> None:
    source = state.u
    cases = (
        (operators.average(source), axis_coefficient_average(source)),
        (operators.primitive(source), axis_coefficient_primitive(source)),
        (operators.parameterPrimitive(source), axis_coefficient_parameter_primitive(source)),
        (operators.mulY(source), axis_coefficient_multiply_y(source)),
        (operators.j1(source), axis_coefficient_regular_inverse(source, 1)),
        (operators.j2(source), axis_coefficient_regular_inverse(source, 2)),
    )
    for actual, expected in cases:
        _same_jets(actual, expected)
        assert actual.epsilon == state.epsilon
        assert actual.paper_exact is False
        assert actual.global_axis_norm_certified is False


def test_bilinear_fields_are_exact_assemblies_of_landed_primitives(state, operators) -> None:
    left = state.u
    right = state.phi
    cases = (
        (operators.product(left, right), axis_coefficient_product(left, right)),
        (operators.param1(left, right), axis_coefficient_inverse_param_product(left, right, 1)),
        (operators.param2(left, right), axis_coefficient_inverse_param_product(left, right, 2)),
        (operators.dot1(left, right), axis_coefficient_inverse_dot_product(left, right, 1)),
        (operators.dot2(left, right), axis_coefficient_inverse_dot_product(left, right, 2)),
        (operators.mixed1(left, right), axis_coefficient_inverse_mixed(left, right, 1)),
        (operators.mixed2(left, right), axis_coefficient_inverse_mixed(left, right, 2)),
    )
    for actual, expected in cases:
        _same_jets(actual, expected, eta=0.21)
        assert actual.epsilon == state.epsilon
        assert actual.paper_exact is False
        assert actual.global_axis_norm_certified is False


def test_snake_case_aliases_do_not_define_new_operators(state, operators) -> None:
    _same_jets(operators.parameter_primitive(state.phi), operators.parameterPrimitive(state.phi))
    _same_jets(operators.mul_y(state.phi), operators.mulY(state.phi))


def test_record_is_anchored_to_actual_reference_epsilon(state, operators) -> None:
    assert operators.epsilon == state.epsilon
    assert operators.reference_origin == "actual SchedulePressure referencePair"
    assert operators.paper_exact is False
    assert operators.global_axis_norm_certified is False

    alien = AxisCoefficientJetState(
        epsilon=2.0 * state.epsilon,
        origin="alien coefficient scale",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        operators.average(alien)
    with pytest.raises(ValueError, match="matching epsilon"):
        operators.product(state.phi, alien)


def test_factory_rejects_non_reference_objects() -> None:
    with pytest.raises(TypeError, match="ActualScheduleReferenceAxisState"):
        actual_schedule_coefficient_operators(object())


def test_operator_record_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_operators.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "coefficientOperators" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
