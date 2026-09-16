from fractions import Fraction as F
import pytest

from openai_ns_reconstruction.kokuno_base_heat_eta_endpoint import (
    heat_eta_endpoint_jet,
    heat_origin_derivative,
    verify_endpoint_pair,
)


def test_exact_fixture_endpoint_jets():
    h = F(1, 200)
    X = F(7, 3)
    minus, plus = verify_endpoint_pair(h, X)
    assert plus.Z == 0
    assert plus.Z_eta == F(-12, 7)
    assert plus.Z_etaeta == F(-12, 7)
    assert plus.H == 1
    assert plus.H_eta == F(603, 70000)
    assert minus.H_eta == -F(603, 70000)
    assert plus.H_etaeta == F(188017209, 4900000000)
    assert minus.H_etaeta == plus.H_etaeta


def test_origin_jet_matches_public_rising_factorial_formula():
    h = F(1, 200)
    assert heat_origin_derivative(h, 0) == 1
    assert heat_origin_derivative(h, 1) == -h * (1 + h)
    assert heat_origin_derivative(h, 2) == h * (1 + h) ** 2 * (2 + h)


def test_endpoint_first_derivative_has_expected_closed_form():
    h = F(1, 211)
    X = F(19, 7)
    plus = heat_eta_endpoint_jet(h, X, 1)
    assert plus.H_eta == 4 * h * (1 + h) / X


def test_endpoint_second_derivative_has_expected_closed_form():
    h = F(1, 211)
    X = F(19, 7)
    plus = heat_eta_endpoint_jet(h, X, 1)
    expected = 16 * h * (1 + h) ** 2 * (2 + h) / X**2 + 4 * h * (1 + h) / X
    assert plus.H_etaeta == expected


def test_transcription_mutation_is_exactly_visible():
    h = F(1, 200)
    X = F(7, 3)
    plus = heat_eta_endpoint_jet(h, X, 1)
    mutated_h1 = heat_origin_derivative(h, 1) + F(1, 2**40)
    mutated = mutated_h1 * plus.Z_eta
    assert mutated != plus.H_eta
    assert mutated - plus.H_eta == F(-3, 7 * 2**38)


def test_wrong_second_coordinate_factor_is_detected():
    h = F(1, 200)
    X = F(7, 3)
    plus = heat_eta_endpoint_jet(h, X, 1)
    H1 = heat_origin_derivative(h, 1)
    H2 = heat_origin_derivative(h, 2)
    wrong_z_etaeta = F(-2, 1) / X
    wrong = H2 * plus.Z_eta**2 + H1 * wrong_z_etaeta
    assert wrong != plus.H_etaeta
    assert wrong - plus.H_etaeta == H1 * F(2, 1) / X


def test_fail_closed_domain_and_exact_type_checks():
    with pytest.raises(ValueError):
        heat_eta_endpoint_jet(F(1, 100), F(1), 1)
    with pytest.raises(ValueError):
        heat_eta_endpoint_jet(F(1, 200), F(0), 1)
    with pytest.raises(ValueError):
        heat_eta_endpoint_jet(F(1, 200), F(1), 0)
    with pytest.raises(TypeError):
        heat_eta_endpoint_jet(0.005, F(1), 1)
    with pytest.raises(TypeError):
        heat_eta_endpoint_jet(F(1, 200), 1.0, 1)
