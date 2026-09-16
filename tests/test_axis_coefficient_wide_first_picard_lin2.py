from decimal import Decimal, localcontext
import math
from types import SimpleNamespace

import pytest

from openai_ns_reconstruction.axis_coefficient_data import (
    actual_schedule_axis_coefficient_data,
)
from openai_ns_reconstruction.axis_coefficient_operators import (
    actual_schedule_coefficient_operators,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_lin2 import (
    ActualScheduleWideFirstPicardLin2State,
    MixedScaleFirstPicardLin2CoefficientJet,
    wide_first_picard_lin2_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96
FIELDS = (
    "ordinary_reference",
    "ordinary_inverse_lambda_numerator",
    "ordinary_inverse_lambda_squared_numerator",
    "pressure_linear_inverse_lambda_numerator",
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def lin2_state(x1):
    return wide_first_picard_lin2_state(x1)


def _u_channels(x1, n: int, m: int, eta: float) -> tuple[Decimal, ...]:
    """Return the four actual x1 input channels, including full pressure jets."""

    _, axial = x1.jet_pair(n, m, eta)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return (
            +axial.reference,
            +axial.inverse_lambda_numerator,
            +axial.inverse_lambda_squared_numerator,
            +(
                x1.remainder.axial.pressure_normalized_factor(n, m, eta)
                / Decimal(2)
            ),
        )


def _row_zero_product(left, right, m: int, eta: float) -> Decimal:
    """Independent eta-Leibniz product for radial-zero coefficient fields."""

    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +sum(
            Decimal(math.comb(m, k))
            * Decimal.from_float(left.jet(0, k, eta))
            * Decimal.from_float(right.jet(0, m - k, eta))
            for k in range(m + 1)
        )


def _linear_coefficient(data, m: int, eta: float) -> Decimal:
    """Compute the pinned radial-zero ``A - 4 A eta uStar + d uStarEta`` jet."""

    with localcontext() as ctx:
        ctx.prec = PRECISION
        A = Decimal.from_float(data.A)
        one = Decimal.from_float(data.one.jet(0, m, eta))
        eta_u_star = _row_zero_product(data.eta, data.uStar, m, eta)
        d_u_star_eta = _row_zero_product(data.d, data.uStarEta, m, eta)
        return +(A * one - Decimal(4) * A * eta_u_star + d_u_star_eta)


def _independent_lin2_oracle(x1, n: int, m: int, eta: float) -> tuple[Decimal, ...]:
    """Direct split-channel convolution, without importing the lin2 branch."""

    if n == 0:
        return (Decimal(0),) * len(FIELDS)

    data = actual_schedule_axis_coefficient_data(x1.reference)
    source_row = n - 1
    divisor = Decimal(n * n)
    result = [Decimal(0) for _ in FIELDS]
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for channel in range(len(FIELDS)):
            total = Decimal(0)
            for k in range(m + 1):
                weight = Decimal(math.comb(m, k))
                source = _u_channels(x1, source_row, m - k, eta)[channel]
                total += weight * _linear_coefficient(data, k, eta) * source

                # dot1(wStar,u): the radial-zero left field leaves the
                # Euler factor source_row in the only surviving slot.
                total += (
                    weight
                    * Decimal(source_row)
                    * Decimal.from_float(data.wStar.jet(0, k, eta))
                    * source
                )

                # param1(u,hStar): eta differentiation acts on the left u.
                total += (
                    weight
                    * _u_channels(x1, source_row, k + 1, eta)[channel]
                    * Decimal.from_float(data.hStar.jet(0, m - k, eta))
                )
            result[channel] = +(total / divisor)
    return tuple(result)


def test_lin2_binds_one_genuine_x1_and_preserves_truth_flags(lin2_state, x1) -> None:
    assert isinstance(lin2_state, ActualScheduleWideFirstPicardLin2State)
    assert lin2_state.x1 is x1
    assert lin2_state.Lambda == x1.Lambda
    assert lin2_state.epsilon == x1.epsilon
    assert lin2_state.lin2_x1_materialized is True
    assert lin2_state.axial_lin2_materialized is True
    assert lin2_state.natural_remainder_x1_materialized is False
    assert lin2_state.fixed_point_materialized is False
    assert lin2_state.paper_exact is False


def test_ordinary_reference_matches_landed_raw_operator_composition(lin2_state, x1) -> None:
    """The Lambda^0 term is the raw lin2 expression, before outer inverseL."""

    n, m, eta = 2, 0, -0.11
    data = actual_schedule_axis_coefficient_data(x1.reference)
    operators = actual_schedule_coefficient_operators(x1.reference)
    reference_u = x1.reference.u

    def axial_linear(i: int, q: int, z: float) -> float:
        return (
            data.A * data.one.jet(i, q, z)
            - 4.0 * data.A * operators.product(data.eta, data.uStar).jet(i, q, z)
            + operators.product(data.d, data.uStarEta).jet(i, q, z)
        )

    linear_state = AxisCoefficientJetState(
        epsilon=x1.epsilon,
        origin="test-only pinned axialLinearCoefficient",
        _jet_provider=axial_linear,
    )
    raw = AxisCoefficientJetState(
        epsilon=x1.epsilon,
        origin="test-only raw lin2 sum",
        _jet_provider=lambda i, q, z: (
            operators.j1(operators.product(linear_state, reference_u)).jet(i, q, z)
            + operators.dot1(data.wStar, reference_u).jet(i, q, z)
            + operators.param1(reference_u, data.hStar).jet(i, q, z)
        ),
    )

    actual = lin2_state.jet(n, m, eta)
    expected = raw.jet(n, m, eta)
    assert math.isclose(
        float(actual.ordinary_reference),
        expected,
        rel_tol=3e-13,
        abs_tol=3e-13,
    )

    # NaturalRemainder applies product(inverseL, lin2) later.  This state is
    # explicitly the raw branch, so an outer inverseL must change this sample.
    outer_inverse = operators.product(data.inverseL, raw).jet(n, m, eta)
    assert not math.isclose(
        float(actual.ordinary_reference),
        outer_inverse,
        rel_tol=1e-10,
        abs_tol=1e-10,
    )


def test_all_four_split_channels_match_independent_differential_oracle(lin2_state, x1) -> None:
    n, m, eta = 3, 0, 0.03
    actual = lin2_state.jet(n, m, eta)
    expected = _independent_lin2_oracle(x1, n, m, eta)
    assert isinstance(actual, MixedScaleFirstPicardLin2CoefficientJet)
    for name, value in zip(FIELDS, expected):
        assert getattr(actual, name) == value
    assert actual.Lambda == x1.Lambda
    assert actual.amplitude_log == (
        x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
    )


def test_full_pressure_jets_and_signed_log_keep_amplitude_derivatives(lin2_state, x1) -> None:
    eta = 0.07
    actual = lin2_state.jet(3, 0, eta)
    assert actual.pressure_linear_inverse_lambda_numerator != 0
    amplitude = x1.remainder.axial.wide_pressure.amplitude
    data = actual_schedule_axis_coefficient_data(x1.reference)
    assert amplitude.log_amplitude(eta - 0.02) != actual.amplitude_log
    assert amplitude.log_amplitude(eta + 0.02) != actual.amplitude_log
    assert amplitude.Lambda * Decimal.from_float(
        data.normalizedGradient.jet(0, 0, eta)
    ) != 0

    term = actual.pressure_linear_term_log()
    numerator = actual.pressure_linear_inverse_lambda_numerator
    if numerator == 0:
        assert term.sign == 0
        assert term.log_scale is None
        assert term.log_factor is None
    else:
        with localcontext() as ctx:
            ctx.prec = PRECISION
            assert term.sign == (1 if numerator > 0 else -1)
            assert term.log_scale == Decimal(2) * actual.amplitude_log
            assert term.log_factor == +(abs(numerator).ln() - actual.Lambda.ln())
        try:
            projected = term.to_binary64()
        except ArithmeticError as exc:
            assert "binary64" in str(exc)
        else:
            assert projected != 0.0
            assert math.isfinite(projected)


def test_row_zero_is_exact_zero_but_keeps_scale_metadata(lin2_state, x1) -> None:
    actual = lin2_state.jet(0, 2, 0.07)
    for name in FIELDS:
        assert getattr(actual, name) == 0
    assert actual.Lambda == x1.Lambda
    assert actual.amplitude_log == (
        x1.remainder.axial.wide_pressure.amplitude.log_amplitude(0.07)
    )
    assert actual.pressure_linear_term_log().sign == 0


def test_scale_helpers_keep_two_ordinary_powers_and_one_pressure_term(lin2_state) -> None:
    actual = lin2_state.jet(3, 0, 0.03)
    first, second = actual.ordinary_correction_terms_decimal()
    with localcontext() as ctx:
        ctx.prec = PRECISION
        assert first == +(actual.ordinary_inverse_lambda_numerator / actual.Lambda)
        assert second == +(
            actual.ordinary_inverse_lambda_squared_numerator
            / actual.Lambda
            / actual.Lambda
        )
    assert actual.pressure_linear_terms_log() == (actual.pressure_linear_term_log(),)


def _poly_derivative(coefficients: tuple[Decimal, ...], order: int, eta: Decimal) -> Decimal:
    return sum(
        (
            Decimal(math.factorial(power))
            / Decimal(math.factorial(power - order))
            * coefficient
            * eta ** (power - order)
        )
        for power, coefficient in enumerate(coefficients)
        if power >= order
    )


# Test-only finite polynomial data.  These numbers are structural diagnostics,
# not SchedulePressure data and never stand in for theorem evidence.
_SYNTHETIC_Q = (Decimal(2), Decimal(-1), Decimal(1))
_SYNTHETIC_W = (Decimal(-1), Decimal(2), Decimal(1))
_SYNTHETIC_H = (Decimal(1), Decimal(-1), Decimal(2))
_SYNTHETIC_ORDINARY = (
    (
        (Decimal(1), Decimal(2)),
        (Decimal(2), Decimal(-1), Decimal(1)),
        (Decimal(3), Decimal(1)),
    ),
    (
        (Decimal(-1), Decimal(1)),
        (Decimal(1), Decimal(1), Decimal(1)),
        (Decimal(2), Decimal(-1)),
    ),
    (
        (Decimal(2), Decimal(-1)),
        (Decimal(3), Decimal(1)),
        (Decimal(1), Decimal(2)),
    ),
)
_SYNTHETIC_PRESSURE = (
    (Decimal(1), Decimal(-1), Decimal(1)),
    (Decimal(2), Decimal(1)),
    (Decimal(-1), Decimal(2), Decimal(1)),
)


def _synthetic_pressure_jet(n: int, m: int, eta: Decimal, amplitude_slope: Decimal) -> Decimal:
    """Derivative of ``a(eta)^2 P(eta)`` after dividing by ``a(eta)^2``."""

    return sum(
        Decimal(math.comb(m, q))
        * (Decimal(2) * amplitude_slope) ** q
        * _poly_derivative(_SYNTHETIC_PRESSURE[n], m - q, eta)
        for q in range(m + 1)
    )


class _SyntheticJetState:
    """Test-only radial-zero jet provider for the production linear action."""

    epsilon = 0.05

    def __init__(self, coefficients: tuple[Decimal, ...]):
        self._coefficients = coefficients

    def jet(self, n: int, m: int, eta: float) -> float:
        if n != 0:
            return 0.0
        return float(_poly_derivative(self._coefficients, m, Decimal(str(eta))))


def test_production_linear_action_matches_expanded_synthetic_polynomials(x1) -> None:
    """Exercise production ``_linear_action`` with explicit diagnostic polynomials."""

    state = wide_first_picard_lin2_state(x1)
    synthetic_data = SimpleNamespace(
        wStar=_SyntheticJetState(_SYNTHETIC_W),
        hStar=_SyntheticJetState(_SYNTHETIC_H),
    )
    amplitude_slope = Decimal("1.5")

    def source(component: str, n: int, m: int, eta: float) -> Decimal:
        point = Decimal(str(eta))
        if component == "pressure":
            return _synthetic_pressure_jet(n, m, point, amplitude_slope)
        channel = FIELDS.index(
            "ordinary_reference"
            if component == "reference"
            else "ordinary_inverse_lambda_numerator"
            if component == "inverse_lambda"
            else "ordinary_inverse_lambda_squared_numerator"
        )
        return _poly_derivative(_SYNTHETIC_ORDINARY[channel][n], m, point)

    original_data = state.axis_data
    original_source = state._source_component
    original_linear = state._axial_linear_jet
    object.__setattr__(state, "axis_data", synthetic_data)
    object.__setattr__(
        state,
        "_source_component",
        source,
    )
    object.__setattr__(
        state,
        "_axial_linear_jet",
        lambda m, eta: _poly_derivative(_SYNTHETIC_Q, m, Decimal(str(eta))),
    )
    try:
        actual = tuple(state._linear_action(component, 3, 3, 0.2) for component in (
            "reference",
            "inverse_lambda",
            "inverse_lambda_squared",
            "pressure",
        ))
    finally:
        object.__setattr__(state, "axis_data", original_data)
        object.__setattr__(state, "_source_component", original_source)
        object.__setattr__(state, "_axial_linear_jet", original_linear)

    # These exact Fraction evaluations are expanded independently from the
    # production loop: [(q + 2 w) f + h f'] / 9, differentiated three times
    # at eta = 1/5.  Pressure uses (D + 3)^m P for variable a(eta)^2.
    expected = (
        Decimal(2),
        Decimal(-2),
        Decimal(4),
        Decimal(121616) / Decimal(1875),
    )
    for actual_value, expected_value in zip(actual, expected):
        assert math.isclose(
            float(actual_value),
            float(expected_value),
            rel_tol=2e-14,
            abs_tol=2e-14,
        )


def test_guards_reject_bad_indices_eta_and_surrogate_state(lin2_state) -> None:
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        lin2_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        lin2_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        lin2_state.jet(1, 0, 2.0)
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_lin2_state(object())
