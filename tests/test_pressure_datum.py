import math

import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.natural_axis_range import (
    NaturalAxisRangeParameters,
    locate_unique_H_root,
)
from openai_ns_reconstruction.pressure_datum import IdealPrefixPressureDatum


def test_minimal_ideal_prefix_witness_uses_exact_lean_threshold():
    datum = IdealPrefixPressureDatum.minimal()
    assert datum.B == 2.0
    assert datum.cap == 1.0
    assert datum.total_mass == 20.0
    assert datum.positive_exponent_mass == 20.0
    assert datum.least_negative_pressure_on_unit_interval == -2.5

    with pytest.raises(ValueError):
        IdealPrefixPressureDatum(B=1.999)
    with pytest.raises(ValueError):
        IdealPrefixPressureDatum(B=float("nan"))


def test_weight_is_exact_ideal_prefix_and_zero_continuation_is_integrable():
    datum = IdealPrefixPressureDatum.minimal()
    for y in (-25.0, -1.0, 0.0):
        assert datum.weight(y) == pytest.approx(4.0 * math.exp(y / 5.0))
        assert datum.exponent(y) == 1.0
    assert datum.weight(1e-12) == 0.0
    assert datum.weight(10.0) == 0.0
    assert datum.exponent(10.0) == 1.0

    mass, error = quad(datum.weight, -math.inf, math.inf, epsabs=1e-12, epsrel=1e-12)
    assert error < 1e-9
    assert mass == pytest.approx(datum.total_mass, rel=1e-12, abs=1e-12)


def test_closed_form_pressure_matches_integral_and_derivative_signs():
    datum = IdealPrefixPressureDatum.minimal()
    for eta in (-1.0, -0.25, 0.0, 0.4, 1.0):
        kernel = datum.kernel(eta)
        integral, _ = quad(
            lambda y: datum.weight(y) * kernel,
            -math.inf,
            math.inf,
            epsabs=1e-12,
            epsrel=1e-12,
        )
        assert datum.pressure(eta) == pytest.approx(-0.5 * integral, rel=1e-12, abs=1e-12)
        P, dP = datum.pressure_data_point(eta)
        assert P <= -1.0
        if eta < 0.0:
            assert dP < 0.0
        elif eta > 0.0:
            assert dP > 0.0
        else:
            assert dP == 0.0

    eps = 1e-6
    eta = -0.31
    finite_difference = (datum.pressure(eta + eps) - datum.pressure(eta - eps)) / (2.0 * eps)
    assert finite_difference == pytest.approx(datum.pressure_derivative(eta), rel=2e-10, abs=2e-10)


def test_natural_axis_root_has_the_formal_positive_Z_separation_for_witness():
    params = NaturalAxisRangeParameters(h=0.005, j=0.02)
    root = locate_unique_H_root(params, atol=1e-15)
    datum = IdealPrefixPressureDatum.minimal()

    z_at_root = datum.natural_axis_Z(params.h, params.j, root.eta0)
    assert z_at_root > params.j / 5.0
    assert z_at_root > 2.0 * params.cutoff_delta


def test_pressure_datum_guards_nonfinite_and_out_of_domain_inputs():
    datum = IdealPrefixPressureDatum.minimal()
    with pytest.raises(ValueError):
        datum.weight(float("inf"))
    with pytest.raises(ValueError):
        datum.exponent(float("-inf"))
    with pytest.raises(ValueError):
        datum.kernel(float("nan"))
    with pytest.raises(ValueError):
        datum.pressure_data_point(1.0001)
