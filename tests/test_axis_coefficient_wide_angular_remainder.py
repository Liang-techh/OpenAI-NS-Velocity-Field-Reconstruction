from decimal import Decimal, localcontext
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_amplitude import (
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.axis_coefficient_natural_remainder import (
    actual_schedule_natural_remainder,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import AxisCoefficientJetState
from openai_ns_reconstruction.axis_coefficient_wide_angular_remainder import (
    MixedScaleAngularCoefficientJet,
    actual_schedule_reference_wide_angular_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def angular():
    return actual_schedule_reference_wide_angular_remainder_state(_schedule_data(), 0.05)


def _constant_state(epsilon: float, value: float) -> AxisCoefficientJetState:
    def provider(n: int, m: int, eta: float) -> float:
        del eta
        return value if n == 0 and m == 0 else 0.0

    return AxisCoefficientJetState(
        epsilon=epsilon,
        origin="regression-only exact constant coefficient state",
        _jet_provider=provider,
    )


def test_state_is_bound_to_actual_reference_resolvent_and_lambda(angular) -> None:
    assert angular.epsilon == angular.reference.epsilon
    assert angular.operators.epsilon == angular.epsilon
    assert angular.axis_data.epsilon == angular.epsilon
    assert angular.resolvent.epsilon == angular.epsilon
    assert angular.resolvent.coefficientwise_series_exact is True
    assert angular.Lambda == actual_schedule_amplitude_log_state(
        _schedule_data(), 0.05
    ).Lambda
    assert angular.mixed_scale_angular_remainder_complete is True
    assert angular.paper_exact is False
    assert angular.global_axis_norm_certified is False


def test_ordinary_base_matches_existing_remainder_at_t_zero(angular) -> None:
    # t=0 and a=0 are regression oracles only.  The angular component is
    # independent of a; these values isolate resolvent(inverseL*(lin1+quad1))
    # in the independently landed full naturalRemainder composition.
    zero = _constant_state(angular.epsilon, 0.0)
    remainder = actual_schedule_natural_remainder(angular.reference)
    expected = remainder.apply(
        0.0,
        zero,
        (angular.reference.phi, angular.reference.u),
    )[0]

    eta = -0.17
    for n in range(5):
        for m in range(2):
            assert angular.ordinary_base.jet(n, m, eta) == pytest.approx(
                expected.jet(n, m, eta), rel=3e-11, abs=3e-12
            )


def test_ordinary_slow_matches_t_difference_after_natural_resolvent(angular) -> None:
    # Linearity of inverseL multiplication and the pinned naturalResolvent gives
    # R_phi(0)-R_phi(1)=resolvent(inverseL*slow1).  t=0/1 are regression-only
    # isolating values and are never manuscript scale choices.
    zero = _constant_state(angular.epsilon, 0.0)
    remainder = actual_schedule_natural_remainder(angular.reference)
    pair = (angular.reference.phi, angular.reference.u)
    at_zero = remainder.apply(0.0, zero, pair)[0]
    at_one = remainder.apply(1.0, zero, pair)[0]

    eta = 0.11
    for n in range(5):
        for m in range(2):
            expected = at_zero.jet(n, m, eta) - at_one.jet(n, m, eta)
            assert angular.ordinary_slow.jet(n, m, eta) == pytest.approx(
                expected, rel=5e-10, abs=5e-11
            )


def test_mixed_jet_preserves_nonzero_inverse_lambda_term(angular) -> None:
    eta = 0.0
    location = next(
        (
            (n, m)
            for n in range(12)
            for m in range(2)
            if angular.ordinary_slow.jet(n, m, eta) != 0.0
        ),
        None,
    )
    assert location is not None
    n, m = location

    jet = angular.jet(n, m, eta)
    assert isinstance(jet, MixedScaleAngularCoefficientJet)
    assert jet.Lambda == angular.Lambda
    assert jet.ordinary_base == Decimal.from_float(angular.ordinary_base.jet(n, m, eta))
    assert jet.inverse_lambda_numerator == -Decimal.from_float(
        angular.ordinary_slow.jet(n, m, eta)
    )

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = +(jet.inverse_lambda_numerator / jet.Lambda)
    assert jet.inverse_lambda_term_decimal() == expected
    assert expected != 0


def test_invalid_indices_and_eta_fail_closed(angular) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        angular.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        angular.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        angular.jet(0, 0, 2.0)


def test_wide_angular_remainder_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_angular_remainder.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "lin1" in layer["capability"]
    assert "naturalResolvent" in layer["capability"]
    assert "1/Lambda" in layer["capability"]
    assert "naturalRemainder(x0)" in layer["remaining_boundary"]
    assert "Picard x1" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
