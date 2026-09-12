from decimal import Decimal
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_amplitude import (
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.natural_axis import real_gradient
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.stage1_scale_chain_wide import (
    diagnose_actual_schedule_scale_chain_wide,
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def amplitude():
    return actual_schedule_amplitude_log_state(_schedule_data(), 0.05)


def test_amplitude_is_bound_to_actual_theorem_scale(amplitude) -> None:
    diagnostic = diagnose_actual_schedule_scale_chain_wide(_schedule_data(), 0.05)
    assert diagnostic.upstream_obstruction is None
    assert diagnostic.scale is not None
    assert amplitude.Lambda == diagnostic.scale.Lambda
    assert amplitude.C_exponent == diagnostic.scale.C.exponent_upper
    assert amplitude.epsilon == amplitude.reference.epsilon
    assert amplitude.radial_degree_zero is True
    assert amplitude.paper_exact is False
    assert amplitude.global_axis_norm_certified is False
    assert amplitude.phase_is_numerically_evaluated is True


def test_log_amplitude_at_origin_uses_symbolic_C_without_exponentiating(amplitude) -> None:
    # realPhase(0)=0 exactly in the landed implementation, so this identity
    # checks the theorem-selected normalization without ever materializing C.
    assert amplitude.log_amplitude(0.0) == -amplitude.C_exponent
    jet = amplitude.jet_log(0, 0, 0.0)
    assert jet.sign == 1
    assert jet.log_abs == -amplitude.C_exponent


def test_first_parameter_jet_uses_real_gradient_identity(amplitude) -> None:
    eta = -0.17
    value = amplitude.jet_log(0, 0, eta)
    derivative = amplitude.jet_log(0, 1, eta)
    g = real_gradient(
        amplitude.data.h,
        amplitude.data.j,
        amplitude.data.sigma,
        eta,
    )
    expected_factor = amplitude.Lambda * Decimal.from_float(g)

    assert value.sign == 1
    assert derivative.sign == (1 if expected_factor > 0 else -1)
    assert derivative.log_abs is not None
    assert value.log_abs is not None

    # a' = Lambda*g*a, compared in log magnitude so the actual theorem scale
    # need not fit binary64.
    with pytest.raises(ArithmeticError):
        # The origin is deliberately far below binary64 on the current
        # conservative scale.  The projection must not silently manufacture 0.
        amplitude.binary64_state().jet(0, 0, 0.0)

    assert float(derivative.log_abs - value.log_abs) == pytest.approx(
        float(abs(expected_factor).ln()),
        rel=1e-12,
        abs=1e-12,
    )


def test_radial_rows_above_zero_are_exactly_zero(amplitude) -> None:
    for n in (1, 2, 5):
        for m in (0, 1, 3):
            jet = amplitude.jet_log(n, m, 0.11)
            assert jet.sign == 0
            assert jet.log_abs is None
            assert amplitude.binary64_state().jet(n, m, 0.11) == 0.0


def test_invalid_indices_and_nonfinite_eta_fail_closed(amplitude) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        amplitude.jet_log(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        amplitude.jet_log(0, True, 0.0)
    with pytest.raises(ValueError, match="eta must be finite"):
        amplitude.log_amplitude(float("nan"))


def test_amplitude_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_amplitude.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "realAmplitude" in layer["capability"]
    assert "binary64" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
