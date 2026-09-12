from dataclasses import replace
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_reference_state import (
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_natural_remainder import (
    ActualScheduleReferenceWideNaturalRemainderState,
    actual_schedule_reference_wide_natural_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def remainder_x0():
    return actual_schedule_reference_wide_natural_remainder_state(
        _schedule_data(),
        0.05,
    )


def test_pair_uses_one_actual_schedule_and_one_theorem_scale(remainder_x0) -> None:
    angular = remainder_x0.angular
    axial = remainder_x0.axial

    assert remainder_x0.reference is angular.reference
    assert remainder_x0.epsilon == angular.epsilon == axial.epsilon
    assert remainder_x0.Lambda == angular.Lambda == axial.Lambda
    assert angular.reference.reference.data == axial.reference.reference.data == _schedule_data()
    assert angular.reference.reference.j == axial.reference.reference.j == 0.05
    assert angular.reference.reference.sigma == axial.reference.reference.sigma
    assert remainder_x0.mixed_scale_natural_remainder_x0_complete is True
    assert remainder_x0.picard_x1_materialized is False
    assert remainder_x0.paper_exact is False
    assert remainder_x0.global_axis_norm_certified is False


def test_jet_pair_preserves_both_landed_mixed_scale_representations(remainder_x0) -> None:
    # Assembly must not collapse Decimal 1/Lambda or the axial signed-log
    # pressure.  The returned pair is exactly the two independently landed
    # typed representations at the same coordinate.
    n, m, eta = 2, 1, 0.07
    phi_jet, u_jet = remainder_x0.jet_pair(n, m, eta)

    assert phi_jet == remainder_x0.angular.jet(n, m, eta)
    assert u_jet == remainder_x0.axial.jet(n, m, eta)
    assert phi_jet.Lambda == u_jet.Lambda == remainder_x0.Lambda

    assert phi_jet.inverse_lambda_term_decimal() == (
        remainder_x0.angular.jet(n, m, eta).inverse_lambda_term_decimal()
    )
    assert u_jet.inverse_lambda_term_decimal() == (
        remainder_x0.axial.jet(n, m, eta).inverse_lambda_term_decimal()
    )
    assert u_jet.pressure == remainder_x0.axial.jet(n, m, eta).pressure


def test_pair_rejects_actual_reference_from_different_schedule_j(remainder_x0) -> None:
    # The incompatible reference is itself constructed through the real
    # SchedulePressure theorem path; this is a fail-closed compatibility test,
    # not a surrogate production input.
    other_reference = actual_schedule_reference_axis_state(_schedule_data(), 0.04)
    mismatched_axial = replace(remainder_x0.axial, reference=other_reference)

    with pytest.raises(ValueError, match="schedule j mismatch"):
        ActualScheduleReferenceWideNaturalRemainderState(
            angular=remainder_x0.angular,
            axial=mismatched_axial,
        )


def test_factory_rejects_non_schedule_data() -> None:
    with pytest.raises(TypeError, match="data must be TailData"):
        actual_schedule_reference_wide_natural_remainder_state(object(), 0.05)


def test_wide_natural_remainder_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_natural_remainder.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "complete mixed-scale naturalRemainder(x0) pair" in layer["capability"]
    assert "1/(2*Lambda)" in layer["remaining_boundary"]
    assert "Picard x1" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
