from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_natural_remainder import (
    actual_schedule_natural_remainder,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import AxisCoefficientJetState
from openai_ns_reconstruction.axis_coefficient_wide_axial_remainder import (
    MixedScaleAxialCoefficientJet,
    actual_schedule_reference_wide_axial_remainder_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def axial():
    return actual_schedule_reference_wide_axial_remainder_state(_schedule_data(), 0.05)


def _constant_state(epsilon: float, value: float) -> AxisCoefficientJetState:
    def provider(n: int, m: int, eta: float) -> float:
        del eta
        return value if n == 0 and m == 0 else 0.0

    return AxisCoefficientJetState(
        epsilon=epsilon,
        origin="regression-only exact constant coefficient state",
        _jet_provider=provider,
    )


def test_state_is_bound_to_actual_reference_pressure_and_lambda(axial) -> None:
    assert axial.epsilon == axial.reference.epsilon
    assert axial.operators.epsilon == axial.epsilon
    assert axial.axis_data.epsilon == axial.epsilon
    assert axial.wide_pressure.epsilon == axial.epsilon
    assert axial.Lambda == axial.wide_pressure.amplitude.Lambda
    assert axial.mixed_scale_axial_remainder_complete is True
    assert axial.paper_exact is False
    assert axial.global_axis_norm_certified is False


def test_ordinary_base_matches_existing_remainder_with_zero_pressure_and_t(axial) -> None:
    # t=0 and a=0 are regression oracles only.  They isolate lin2 through the
    # independently landed full naturalRemainder composition; they are not
    # manuscript parameter selections.
    zero = _constant_state(axial.epsilon, 0.0)
    remainder = actual_schedule_natural_remainder(axial.reference)
    expected = remainder.apply(
        0.0,
        zero,
        (axial.reference.phi, axial.reference.u),
    )[1]

    eta = -0.17
    for n in range(4):
        for m in range(2):
            assert axial.ordinary_base.jet(n, m, eta) == pytest.approx(
                expected.jet(n, m, eta), rel=2e-12, abs=2e-13
            )


def test_ordinary_slow_matches_t_difference_in_existing_remainder(axial) -> None:
    # With a=0, R_u(t=0)-R_u(t=1) is exactly inverseL*slow2.  Again, t=0/1
    # exist only as an independent decomposition oracle.
    zero = _constant_state(axial.epsilon, 0.0)
    remainder = actual_schedule_natural_remainder(axial.reference)
    pair = (axial.reference.phi, axial.reference.u)
    at_zero = remainder.apply(0.0, zero, pair)[1]
    at_one = remainder.apply(1.0, zero, pair)[1]

    eta = 0.11
    for n in range(4):
        for m in range(2):
            expected = at_zero.jet(n, m, eta) - at_one.jet(n, m, eta)
            assert axial.ordinary_slow.jet(n, m, eta) == pytest.approx(
                expected, rel=3e-11, abs=3e-12
            )


def test_pressure_after_inverse_l_uses_radial_zero_leibniz_identity(axial) -> None:
    # inverseL is one of the landed AxisData fields and is exactly radial-degree
    # zero.  Therefore the final product with the wide pressure has a simple
    # independent row identity.  This avoids the invalid constant-amplitude
    # oracle: parameterPrimitive differentiates the *actual* amplitude, so
    # replacing a(eta) by constant 1 would erase a genuine Lambda-sized term.
    eta = 0.07
    inverse_l = axial.axis_data.inverseL

    for radial_row in range(1, 4):
        for derivative_order in range(3):
            assert inverse_l.jet(radial_row, derivative_order, eta) == 0.0

    for n in range(6):
        with localcontext() as ctx:
            ctx.prec = PRECISION
            l0 = Decimal.from_float(inverse_l.jet(0, 0, eta))
            l1 = Decimal.from_float(inverse_l.jet(0, 1, eta))
            p0 = axial.wide_pressure.normalized_factor(n, 0, eta)
            p1 = axial.wide_pressure.normalized_factor(n, 1, eta)
            expected0 = +(l0 * p0)
            expected1 = +(l0 * p1 + l1 * p0)
        assert axial.pressure_normalized_factor(n, 0, eta) == expected0
        assert axial.pressure_normalized_factor(n, 1, eta) == expected1


def test_mixed_jet_keeps_inverse_lambda_and_pressure_as_separate_scales(axial) -> None:
    eta = 0.0
    nonzero_pressure_row = next(
        (
            n
            for n in range(12)
            if axial.pressure_normalized_factor(n, 0, eta) != 0
        ),
        None,
    )
    assert nonzero_pressure_row is not None

    jet = axial.jet(nonzero_pressure_row, 0, eta)
    assert isinstance(jet, MixedScaleAxialCoefficientJet)
    assert jet.Lambda == axial.Lambda
    assert jet.ordinary_base == Decimal.from_float(
        axial.ordinary_base.jet(nonzero_pressure_row, 0, eta)
    )
    assert jet.inverse_lambda_numerator == -Decimal.from_float(
        axial.ordinary_slow.jet(nonzero_pressure_row, 0, eta)
    )
    assert jet.pressure.sign != 0
    assert jet.pressure.log_abs is not None

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected_slow = +(jet.inverse_lambda_numerator / jet.Lambda)
    assert jet.inverse_lambda_term_decimal() == expected_slow

    # A nonzero pressure term may be too small for binary64, but it may never be
    # silently rounded to zero by this bridge.
    try:
        projected = jet.pressure.to_binary64()
    except ArithmeticError as exc:
        assert "binary64" in str(exc)
    else:
        assert projected != 0.0
        assert math.isfinite(projected)


def test_invalid_indices_and_eta_fail_closed(axial) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        axial.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        axial.pressure_normalized_factor(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        axial.pressure_jet_log(0, 0, 2.0)


def test_wide_axial_remainder_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_axial_remainder.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "lin2" in layer["capability"]
    assert "1/Lambda" in layer["capability"]
    assert "angular" in layer["remaining_boundary"]
    assert "Picard x1" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
