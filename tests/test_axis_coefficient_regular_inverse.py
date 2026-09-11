import json
import math
from pathlib import Path

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.axis_coefficient_reference_state import (
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.axis_coefficient_regular_inverse import (
    axis_coefficient_regular_inverse,
    radial_divisor,
    regular_inverse_jet,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


ROOT = Path(__file__).resolve().parents[1]


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def state():
    return actual_schedule_reference_axis_state(_schedule_data(), 0.05)


def _finite_radial_polynomial(state, degree: int, eta: float):
    coefficients = [state.coefficient(n, eta) for n in range(degree + 1)]

    def value(y: float) -> float:
        return sum(coefficient * y**n for n, coefficient in enumerate(coefficients))

    return value


def _independent_zero_datum_inverse(source, r: int, y: float) -> tuple[float, float]:
    # Independent Green-kernel evaluation of the regular solution of
    # y*g'' + r*g' = source with g(0)=0.  Production uses only the pinned
    # coefficient row map and never this quadrature representation.
    if r == 1:
        def integrand(s: float) -> float:
            if s == 0.0:
                return 0.0
            return math.log(y / s) * source(s)
    else:
        def integrand(s: float) -> float:
            return (1.0 - (s / y) ** (r - 1)) * source(s) / float(r - 1)

    return quad(integrand, 0.0, y, epsabs=1e-12, epsrel=1e-12)


@pytest.mark.parametrize("r", [1, 2])
def test_regular_inverse_matches_independent_green_kernel(state, r: int) -> None:
    inverse = axis_coefficient_regular_inverse(state.phi, r)
    eta = 0.27
    y = 0.41
    degree = 6
    source = _finite_radial_polynomial(state.phi, degree, eta)

    oracle, error = _independent_zero_datum_inverse(source, r, y)
    reconstructed = sum(
        inverse.coefficient(n, eta) * y**n for n in range(degree + 2)
    )
    assert error < 1e-9
    assert reconstructed == pytest.approx(oracle, rel=3e-11, abs=3e-11)


@pytest.mark.parametrize("r", [1, 2, 4])
def test_regular_inverse_satisfies_physical_radial_equation(state, r: int) -> None:
    inverse = axis_coefficient_regular_inverse(state.phi, r)
    eta = -0.31
    y = 0.37
    degree = 6

    source_value = sum(
        state.phi.coefficient(n, eta) * y**n for n in range(degree + 1)
    )
    first = sum(
        n * inverse.coefficient(n, eta) * y ** (n - 1)
        for n in range(1, degree + 2)
    )
    second = sum(
        n * (n - 1) * inverse.coefficient(n, eta) * y ** (n - 2)
        for n in range(2, degree + 2)
    )
    assert y * second + float(r) * first == pytest.approx(
        source_value, rel=3e-12, abs=3e-12
    )


def test_regular_inverse_parameter_jet_matches_finite_difference(state) -> None:
    inverse = axis_coefficient_regular_inverse(state.phi, 2)
    eta = -0.23
    n = 3
    step = 2.0e-5

    def direct_inverse_coefficient(x: float) -> float:
        return state.phi.coefficient(n - 1, x) / radial_divisor(2, n - 1)

    oracle = (
        direct_inverse_coefficient(eta + step)
        - direct_inverse_coefficient(eta - step)
    ) / (2.0 * step)
    assert inverse.jet(n, 1, eta) == pytest.approx(oracle, rel=2e-5, abs=2e-7)


@pytest.mark.parametrize("r, factor", [(1, 0.25), (2, 1.0 / 6.0)])
def test_axial_reference_degree_one_moves_to_degree_two(state, r: int, factor: float) -> None:
    inverse = axis_coefficient_regular_inverse(state.u, r)
    eta = 0.18
    for m in range(4):
        assert inverse.jet(0, m, eta) == 0.0
        assert inverse.jet(1, m, eta) == 0.0
        assert inverse.jet(2, m, eta) == pytest.approx(
            factor * state.u.jet(1, m, eta), rel=2e-14, abs=1e-12
        )
        assert inverse.jet(3, m, eta) == 0.0


def test_regular_inverse_jet_matches_pinned_row_formula(state) -> None:
    eta = -0.39
    for r in (1, 2, 5):
        for m in range(3):
            assert regular_inverse_jet(state.phi, r, 0, m, eta) == 0.0
            for n in range(1, 7):
                assert regular_inverse_jet(state.phi, r, n, m, eta) == pytest.approx(
                    state.phi.jet(n - 1, m, eta) / radial_divisor(r, n - 1),
                    rel=2e-15,
                    abs=1e-14,
                )


def test_regular_inverse_stays_fail_closed_for_axis_space_membership(state) -> None:
    inverse = axis_coefficient_regular_inverse(state.phi, 2)
    assert inverse.paper_exact is False
    assert inverse.global_axis_norm_certified is False
    assert inverse.epsilon == state.epsilon
    assert "AxisOperators regularInverse(r=2)" in inverse.origin


def test_regular_inverse_preserves_guards(state) -> None:
    for bad_r in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            axis_coefficient_regular_inverse(state.phi, bad_r)

    inverse = axis_coefficient_regular_inverse(state.phi, 2)
    with pytest.raises(ValueError, match="nonnegative integer"):
        inverse.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        inverse.jet(True, 0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        inverse.jet(0, -1, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        inverse.jet(0, 0, 1.100001)


def test_regular_inverse_provenance_is_machine_readable_and_fail_closed() -> None:
    manifest = json.loads(
        (
            ROOT
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_regular_inverse.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    assert "AxisOperators.regularInverse" in layer["capability"]
    assert "naturalRemainder" in layer["remaining_boundary"]
    assert "global all-index weighted norm" in layer["remaining_boundary"]
    assert (ROOT / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (ROOT / artifact).is_file()
