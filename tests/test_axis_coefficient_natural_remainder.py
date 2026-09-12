import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_data import actual_schedule_axis_coefficient_data
from openai_ns_reconstruction.axis_coefficient_natural_remainder import (
    actual_schedule_natural_remainder,
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
def remainder(reference):
    return actual_schedule_natural_remainder(reference)


def _scaled(state: AxisCoefficientJetState, scalar: float, label: str) -> AxisCoefficientJetState:
    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=label,
        _jet_provider=lambda n, m, eta: scalar * state.jet(n, m, eta),
    )


def _zero(reference) -> AxisCoefficientJetState:
    return AxisCoefficientJetState(
        epsilon=reference.epsilon,
        origin="test-only zero oracle; not a manuscript amplitude",
        _jet_provider=lambda n, m, eta: 0.0,
    )


def test_factory_is_anchored_to_actual_schedule_chain(reference, remainder) -> None:
    assert remainder.epsilon == reference.epsilon
    assert remainder.operators.epsilon == reference.epsilon
    assert remainder.data.epsilon == reference.epsilon
    assert remainder.resolvent.epsilon == reference.epsilon
    assert remainder.paper_exact is False
    assert remainder.global_axis_norm_certified is False
    assert remainder.actual_amplitude_materialized is False


def test_remainder_row_zero_vanishes_from_pinned_regular_inverses(reference, remainder) -> None:
    data = actual_schedule_axis_coefficient_data(reference)
    result = remainder(0.0, data.one, (reference.phi, reference.u))
    for eta in (-0.21, 0.17):
        for m in (0, 1):
            assert result[0].jet(0, m, eta) == 0.0
            assert result[1].jet(0, m, eta) == 0.0


def test_pressure_source_is_quadratic_in_amplitude_state(reference, remainder) -> None:
    """The only a-dependence is source=(a*a)*(phi*phi), hence pressure is quadratic."""

    data = actual_schedule_axis_coefficient_data(reference)
    zero = _zero(reference)
    twice = _scaled(data.one, 2.0, "test-only twice AxisData.one")
    x = (reference.phi, reference.u)

    r0 = remainder(0.0, zero, x)
    r1 = remainder(0.0, data.one, x)
    r2 = remainder(0.0, twice, x)

    eta = -0.13
    for n in (2, 3):
        # The angular remainder is independent of the pressure amplitude a.
        assert r2[0].jet(n, 0, eta) == pytest.approx(r1[0].jet(n, 0, eta), rel=1e-12, abs=1e-12)
        # Axial pressure contribution scales like a^2.
        delta1 = r1[1].jet(n, 0, eta) - r0[1].jet(n, 0, eta)
        delta2 = r2[1].jet(n, 0, eta) - r0[1].jet(n, 0, eta)
        assert delta2 == pytest.approx(4.0 * delta1, rel=2e-9, abs=2e-10)


def test_remainder_is_affine_in_inverse_lambda_parameter(reference, remainder) -> None:
    """Pinned t-dependence occurs only as -t*slow1 and -t*slow2."""

    data = actual_schedule_axis_coefficient_data(reference)
    x = (reference.phi, reference.u)
    r0 = remainder(0.0, data.one, x)
    r1 = remainder(0.125, data.one, x)
    r2 = remainder(0.25, data.one, x)

    eta = 0.11
    for component in (0, 1):
        for n in (1, 2):
            delta1 = r1[component].jet(n, 0, eta) - r0[component].jet(n, 0, eta)
            delta2 = r2[component].jet(n, 0, eta) - r0[component].jet(n, 0, eta)
            assert delta2 == pytest.approx(2.0 * delta1, rel=2e-9, abs=2e-10)


def test_remainder_fails_closed_on_wrong_scale(reference, remainder) -> None:
    alien = AxisCoefficientJetState(
        epsilon=2.0 * reference.epsilon,
        origin="wrong-scale naturalRemainder regression state",
        _jet_provider=lambda n, m, eta: 0.0,
    )
    with pytest.raises(ValueError, match="matching epsilon"):
        remainder(0.1, alien, (reference.phi, reference.u))
    with pytest.raises(ValueError, match="matching epsilon"):
        remainder(0.1, reference.phi, (reference.phi, alien))


def test_natural_remainder_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_natural_remainder.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "naturalRemainder" in layer["capability"]
    assert "naturalRemainder(x0)" in layer["remaining_boundary"]
    assert "amplitude" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
