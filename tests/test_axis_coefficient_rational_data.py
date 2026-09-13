from decimal import Decimal, localcontext
from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.axis_coefficient_data import (
    actual_schedule_axis_coefficient_data,
)
from openai_ns_reconstruction.axis_coefficient_rational_data import (
    RationalAxisCoefficientData,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def rational_data() -> RationalAxisCoefficientData:
    reference = actual_schedule_reference_axis_state(_data(), 0.05)
    actual = actual_schedule_axis_coefficient_data(reference)
    return RationalAxisCoefficientData(actual)


def _parameters(
    state: RationalAxisCoefficientData,
    eta: float,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    e = Fraction.from_float(eta)
    d = Fraction(1) - e * e
    U = Fraction(4) * e + state.j
    H = state.D * e + d * U
    L = Fraction(1) - Fraction(2) * state.h * e * e
    return e, d, U, H, L, state.sigma


def test_exact_scalar_and_reference_formulas(rational_data) -> None:
    eta = 0.23
    e, _, _, H, L, sigma = _parameters(rational_data, eta)
    assert rational_data.A == Fraction(1, 2) + rational_data.h
    assert rational_data.D == Fraction(1, 2) - rational_data.h

    H_derivatives = (
        H,
        rational_data.D + 4 - 2 * rational_data.j * e - 12 * e * e,
        -2 * rational_data.j - 24 * e,
        Fraction(-24),
    )
    for order, expected in enumerate(H_derivatives):
        assert rational_data.jet_fraction("hStar", 0, order, eta) == expected

    inverse_derivatives = (
        Fraction(1, 1) / L,
        Fraction(4) * rational_data.h * e / L**2,
        Fraction(4) * rational_data.h / L**2
        + Fraction(32) * rational_data.h**2 * e**2 / L**3,
        Fraction(96) * rational_data.h**2 * e / L**3
        + Fraction(384) * rational_data.h**3 * e**3 / L**4,
    )
    for order, expected in enumerate(inverse_derivatives):
        assert rational_data.jet_fraction("inverseL", 0, order, eta) == expected

    chi = H * H / (H * H + sigma * sigma)
    gradient = -(L * H) / (H * H + sigma * sigma)
    assert rational_data.jet_fraction("chi", 0, 0, eta) == chi
    assert rational_data.jet_fraction("normalizedGradient", 0, 0, eta) == gradient
    assert rational_data.angular_reference_jet_fraction(1, 0, eta) == -chi / 4
    assert rational_data.angular_reference_jet_fraction(2, 0, eta) == chi * chi / 48

    exact = rational_data.jet_fraction("hStar", 0, 2, eta)
    with localcontext() as context:
        context.prec = 96
        expected_decimal = +(Decimal(exact.numerator) / Decimal(exact.denominator))
    assert rational_data.jet_decimal("hStar", 0, 2, eta) == expected_decimal


def test_directed_enclosures_are_exact_and_context_independent(rational_data) -> None:
    eta = 0.23
    requests = (
        ("hStar", 12, 0),
        ("inverseL", 12, 0),
        ("chi", 12, 0),
        ("normalizedGradient", 12, 0),
    )
    for name, order, radial in requests:
        exact = rational_data.jet_fraction(name, radial, order, eta)
        with localcontext() as context:
            context.prec = 8
            low_precision = rational_data.jet_enclosure(name, radial, order, eta)
        with localcontext() as context:
            context.prec = 120
            high_precision = rational_data.jet_enclosure(name, radial, order, eta)
        assert Fraction(low_precision.lower) <= exact <= Fraction(low_precision.upper)
        assert low_precision.lower == high_precision.lower
        assert low_precision.upper == high_precision.upper

    exact_reference = rational_data.angular_reference_jet_fraction(2, 12, eta)
    enclosure = rational_data.angular_reference_jet_enclosure(2, 12, eta)
    assert Fraction(enclosure.lower) <= exact_reference <= Fraction(enclosure.upper)


def test_radial_and_domain_guards_fail_closed(rational_data) -> None:
    eta = 0.23
    fixed_names = (
        "one",
        "eta",
        "d",
        "inverseL",
        "uStar",
        "uStarEta",
        "wStar",
        "hStar",
        "normalizedGradient",
        "chi",
    )
    for name in fixed_names:
        assert rational_data.jet_fraction(name, 1, 0, eta) == Fraction(0)
    with pytest.raises(NotImplementedError, match="zStar"):
        rational_data.jet_fraction("zStar", 0, 0, eta)

    inside_endpoint = math.nextafter(1.1, 0.0)
    assert rational_data.jet_fraction("one", 0, 0, inside_endpoint) == Fraction(1)
    with pytest.raises(ValueError, match="pinned window"):
        rational_data.jet_fraction("one", 0, 0, 1.1)
    assert rational_data.paper_exact is False
    assert rational_data.global_axis_norm_certified is False
