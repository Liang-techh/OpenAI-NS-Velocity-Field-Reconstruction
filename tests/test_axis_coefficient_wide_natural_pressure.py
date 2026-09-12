from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_natural_pressure import (
    actual_schedule_reference_wide_pressure_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def pressure():
    return actual_schedule_reference_wide_pressure_log_state(_schedule_data(), 0.05)


def _d(value: float) -> Decimal:
    return Decimal.from_float(float(value))


def _ordinary_times_wide(ordinary, wide, n: int, m: int, eta: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        total = Decimal(0)
        for i in range(n + 1):
            j = n - i
            for k in range(m + 1):
                total += (
                    Decimal(math.comb(m, k))
                    * wide(i, k, eta)
                    * _d(ordinary.jet(j, m - k, eta))
                )
        return +total


def test_pressure_is_bound_to_actual_source_and_scale(pressure) -> None:
    assert pressure.epsilon == pressure.source.epsilon
    assert pressure.amplitude is pressure.source.amplitude
    assert pressure.data is pressure.source.amplitude.data
    assert pressure.pressure_chain_complete is True
    assert pressure.paper_exact is False
    assert pressure.global_axis_norm_certified is False


def test_linear_radial_maps_preserve_the_common_amplitude_square(pressure) -> None:
    eta = -0.13
    for m in range(2):
        assert pressure.primitive_factor(0, m, eta) == 0
        assert pressure.parameter_primitive_factor(0, m, eta) == 0
        assert pressure.multiply_y_factor(0, m, eta) == 0
        for n in range(1, 4):
            with localcontext() as ctx:
                ctx.prec = PRECISION
                source_nm = pressure.source.normalized_factor(n - 1, m, eta)
                source_next = pressure.source.normalized_factor(n - 1, m + 1, eta)
                assert pressure.primitive_factor(n, m, eta) == +(source_nm / Decimal(n))
                assert pressure.parameter_primitive_factor(n, m, eta) == +(
                    source_next / Decimal(n)
                )
                assert pressure.multiply_y_factor(n, m, eta) == source_nm


def test_pressure_input_matches_literal_pinned_chain(pressure) -> None:
    eta = 0.07
    d = pressure.data
    for n in range(4):
        for m in range(2):
            eta_primitive = _ordinary_times_wide(
                d.eta, pressure.primitive_factor, n, m, eta
            )
            d_parameter_primitive = _ordinary_times_wide(
                d.d, pressure.parameter_primitive_factor, n, m, eta
            )
            eta_multiply_y = _ordinary_times_wide(
                d.eta, pressure.multiply_y_factor, n, m, eta
            )
            with localcontext() as ctx:
                ctx.prec = PRECISION
                expected = +(
                    -Decimal(4) * _d(d.A) * eta_primitive
                    + d_parameter_primitive
                    - Decimal(2) * eta_multiply_y
                )
            assert pressure.pressure_input_factor(n, m, eta) == expected


def test_final_j1_uses_exact_regular_inverse_row_divisor(pressure) -> None:
    eta = -0.19
    for m in range(2):
        assert pressure.normalized_factor(0, m, eta) == 0
        for n in range(1, 5):
            with localcontext() as ctx:
                ctx.prec = PRECISION
                expected = +(
                    pressure.pressure_input_factor(n - 1, m, eta)
                    / Decimal(n * n)
                )
            assert pressure.normalized_factor(n, m, eta) == expected


def test_nonzero_pressure_jet_never_silently_becomes_binary64_zero(pressure) -> None:
    eta = 0.0
    nonzero = next(
        (n for n in range(10) if pressure.normalized_factor(n, 0, eta) != 0),
        None,
    )
    assert nonzero is not None
    jet = pressure.jet_log(nonzero, 0, eta)
    assert jet.sign != 0
    assert jet.log_abs is not None
    try:
        projected = pressure.binary64_state().jet(nonzero, 0, eta)
    except ArithmeticError as exc:
        assert "binary64" in str(exc)
    else:
        assert projected != 0.0
        assert math.isfinite(projected)


def test_invalid_indices_eta_and_epsilon_mismatch_fail_closed(pressure) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        pressure.normalized_factor(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        pressure.pressure_input_factor(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        pressure.jet_log(0, 0, 2.0)


def test_wide_pressure_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_natural_pressure.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "pressure" in layer["capability"]
    assert "naturalRemainder(x0)" in layer["remaining_boundary"]
    assert "t=1/Lambda" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
