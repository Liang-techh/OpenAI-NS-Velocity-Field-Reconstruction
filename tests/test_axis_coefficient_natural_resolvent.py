import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_data import actual_schedule_axis_coefficient_data
from openai_ns_reconstruction.axis_coefficient_natural_resolvent import (
    actual_schedule_natural_resolvent,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


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
def resolvent(reference):
    return actual_schedule_natural_resolvent(reference)


def test_natural_resolvent_exact_term_count_is_radial_filtration() -> None:
    assert resolvent_term_counts() == (1, 2, 3, 4, 5, 6)


def resolvent_term_counts() -> tuple[int, ...]:
    from openai_ns_reconstruction.axis_coefficient_natural_resolvent import (
        AxisCoefficientNaturalResolvent,
    )

    return tuple(AxisCoefficientNaturalResolvent.exact_term_count(n) for n in range(6))


def test_natural_resolvent_of_one_recovers_actual_reference_phi(reference, resolvent) -> None:
    """Pinned AxisReference: naturalResolvent(one) is the actual angular reference."""

    data = actual_schedule_axis_coefficient_data(reference)
    resolved = resolvent(data.one)
    for eta in (-0.29, -0.013, 0.24):
        for n in range(7):
            for m in (0, 1, 2):
                assert resolved.jet(n, m, eta) == pytest.approx(
                    reference.phi.jet(n, m, eta),
                    rel=2e-9,
                    abs=3e-13,
                )


def test_natural_resolvent_satisfies_one_plus_q_equation_on_actual_axial_state(
    reference,
    resolvent,
) -> None:
    """Check R(A) + Q(R(A)) = A on a nontrivial actual coefficient state."""

    resolved = resolvent(reference.u)
    q_resolved = resolvent.operator(resolved)
    for eta in (-0.21, 0.17):
        for n in range(6):
            for m in (0, 1):
                lhs = resolved.jet(n, m, eta) + q_resolved.jet(n, m, eta)
                rhs = reference.u.jet(n, m, eta)
                assert lhs == pytest.approx(rhs, rel=2e-9, abs=5e-13)


def test_natural_resolvent_row_zero_is_exact_identity(reference, resolvent) -> None:
    resolved = resolvent(reference.u)
    for eta in (-0.3, 0.0, 0.31):
        for m in range(3):
            assert resolved.jet(0, m, eta) == reference.u.jet(0, m, eta)


def test_natural_resolvent_fails_closed_on_wrong_scale(reference, resolvent) -> None:
    wrong = AxisCoefficientJetState(
        epsilon=2.0 * reference.epsilon,
        origin="wrong-scale naturalResolvent regression state",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        resolvent(wrong)


def test_natural_resolvent_truth_boundary_is_explicit(reference, resolvent) -> None:
    assert resolvent.epsilon == reference.epsilon
    assert resolvent.paper_exact is False
    assert resolvent.global_axis_norm_certified is False
    assert resolvent.coefficientwise_series_exact is True
    assert resolvent.natural_remainder_materialized is False


def test_natural_resolvent_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_natural_resolvent.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "naturalResolvent" in layer["capability"]
    assert "coefficientwise" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
