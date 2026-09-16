from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.axis_coefficient_wide_natural_source import (
    actual_schedule_reference_natural_source_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def source():
    return actual_schedule_reference_natural_source_log_state(_schedule_data(), 0.05)


def _d(value: float) -> Decimal:
    return Decimal.from_float(float(value))


def _nearest96(value: Fraction) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = ROUND_HALF_EVEN
        return +(Decimal(value.numerator) / Decimal(value.denominator))


def _a_square_bell_fraction(source, order: int, eta: float) -> Fraction:
    """Independent low-order Bell factors for ``a^2``."""

    if order == 0:
        return Fraction(1)
    rational = source.amplitude.rational_data
    t = 2 * Fraction(source.amplitude.Lambda)
    g0 = rational.jet_fraction("normalizedGradient", 0, 0, eta)
    if order == 1:
        return t * g0
    assert order == 2
    g1 = rational.jet_fraction("normalizedGradient", 0, 1, eta)
    return t * g1 + t * t * g0 * g0


def test_source_is_bound_to_actual_reference_pair_and_theorem_amplitude(source) -> None:
    assert source.epsilon == source.amplitude.reference.epsilon
    assert source.angular_square.epsilon == source.epsilon
    assert source.x_is_reference_pair is True
    assert source.radial_amplitude_degree_zero is True
    assert source.paper_exact is False
    assert source.global_axis_norm_certified is False
    assert "productFamily" in source.angular_square.origin or "product" in source.angular_square.origin


def test_zeroth_eta_jet_extracts_exact_common_amplitude_square_scale(source) -> None:
    eta = -0.13
    for n in range(4):
        expected_factor = _d(source.angular_square.jet(n, 0, eta))
        assert source.normalized_factor(n, 0, eta) == expected_factor
        jet = source.jet_log(n, 0, eta)
        if expected_factor == 0:
            assert jet.sign == 0
            assert jet.log_abs is None
            continue
        assert jet.sign == (1 if expected_factor > 0 else -1)
        with localcontext() as ctx:
            ctx.prec = PRECISION
            expected_scale = +(Decimal(2) * source.amplitude.log_amplitude(eta))
            expected_log_factor = +abs(expected_factor).ln()
        assert jet.log_scale == expected_scale
        assert jet.log_factor == expected_log_factor


def test_first_eta_jet_matches_direct_product_rule(source) -> None:
    eta = -0.17
    n = 0
    payload0 = _d(source.angular_square.jet(n, 0, eta))
    payload1 = _d(source.angular_square.jet(n, 1, eta))
    bell0 = _nearest96(_a_square_bell_fraction(source, 0, eta))
    bell1 = _nearest96(_a_square_bell_fraction(source, 1, eta))

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = Decimal(0)
        expected += bell0 * payload1
        expected += bell1 * payload0
        expected = +expected

    assert source.normalized_factor(n, 1, eta) == expected


def test_second_eta_jet_matches_explicit_bell_leibniz_identity(source) -> None:
    eta = 0.09
    n = 1
    p0 = _d(source.angular_square.jet(n, 0, eta))
    p1 = _d(source.angular_square.jet(n, 1, eta))
    p2 = _d(source.angular_square.jet(n, 2, eta))
    bell0 = _nearest96(_a_square_bell_fraction(source, 0, eta))
    bell1 = _nearest96(_a_square_bell_fraction(source, 1, eta))
    bell2 = _nearest96(_a_square_bell_fraction(source, 2, eta))

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = Decimal(0)
        expected += bell0 * p2
        expected += Decimal(2) * bell1 * p1
        expected += bell2 * p0
        expected = +expected

    assert source.normalized_factor(n, 2, eta) == expected


def test_current_nonzero_source_refuses_binary64_underflow(source) -> None:
    eta = 0.0
    nonzero_row = next(
        n for n in range(8) if source.normalized_factor(n, 0, eta) != 0
    )
    jet = source.jet_log(nonzero_row, 0, eta)
    assert jet.sign != 0
    assert jet.log_abs is not None
    assert jet.log_abs < Decimal.from_float(math.ulp(0.0)).ln()
    with pytest.raises(ArithmeticError, match="underflows binary64"):
        source.binary64_state().jet(nonzero_row, 0, eta)


def test_invalid_indices_and_eta_fail_closed(source) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        source.normalized_factor(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        source.jet_log(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        source.jet_log(0, 0, 2.0)


def test_wide_source_provenance_is_machine_readable_and_fail_closed() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_natural_source.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "a^2*phi0^2" in layer["capability"]
    assert "naturalRemainder(x0)" in layer["remaining_boundary"]
    assert "t=1/Lambda" in layer["remaining_boundary"]
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()
