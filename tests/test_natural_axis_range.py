import math

import pytest

from openai_ns_reconstruction.natural_axis import H, chi
from openai_ns_reconstruction.natural_axis_range import (
    NaturalAxisRangeParameters,
    cutoff_parameters_from_margin,
    locate_unique_H_root,
    verify_chi_margin_point,
)


def test_full_printed_parameter_range_and_exact_delta_choice():
    edge = NaturalAxisRangeParameters(h=1.0 / 100.0, j=1.0 / 20.0)
    assert edge.root_bracket == (-edge.j / 4.0, -edge.j / 5.0)
    assert edge.cutoff_delta == edge.j / 10.0

    for h, j in [(0.0, 0.01), (0.011, 0.01), (0.001, 0.0), (0.001, 0.051)]:
        with pytest.raises(ValueError):
            NaturalAxisRangeParameters(h=h, j=j)


def test_formal_root_bracket_has_sign_change_across_range_samples():
    for h, j in [(1e-6, 1e-6), (0.001, 0.001), (0.005, 0.02), (0.01, 0.05)]:
        p = NaturalAxisRangeParameters(h=h, j=j)
        left, right = p.root_bracket
        assert H(h, j, left) < 0.0
        assert H(h, j, right) > 0.0


def test_bisection_materializes_numerical_location_of_unique_H_root():
    p = NaturalAxisRangeParameters(h=0.005, j=0.02)
    witness = locate_unique_H_root(p, atol=1e-15)
    left, right = witness.bracket
    assert left < witness.eta0 < right
    assert abs(witness.residual) <= 1e-14
    assert abs(H(p.h, p.j, witness.eta0)) <= 1e-14


def test_sigma_from_certified_margin_matches_lean_algebra():
    p = NaturalAxisRangeParameters(h=0.005, j=0.02)
    # Use a genuine pointwise H^2 value only to test the algebra.  The module
    # does not claim this establishes the required uniform low-|Z| margin.
    h2 = H(p.h, p.j, 0.0) ** 2
    margin = h2 / 2.0
    witness = cutoff_parameters_from_margin(p, margin)

    assert witness.delta == pytest.approx(p.j / 10.0)
    assert witness.sigma == pytest.approx(math.sqrt(margin) / 20.0)
    assert witness.guaranteed_chi_lower_bound == pytest.approx(400.0 / 401.0)
    assert witness.guaranteed_chi_lower_bound > 0.99

    got = verify_chi_margin_point(p, witness, 0.0)
    assert got == pytest.approx(chi(p.h, p.j, witness.sigma, 0.0))
    assert got >= witness.guaranteed_chi_lower_bound


def test_margin_adapter_fails_closed_without_required_evidence():
    p = NaturalAxisRangeParameters(h=0.005, j=0.02)
    with pytest.raises(ValueError):
        cutoff_parameters_from_margin(p, 0.0)
    with pytest.raises(ValueError):
        cutoff_parameters_from_margin(p, float("nan"))

    witness = cutoff_parameters_from_margin(p, H(p.h, p.j, 0.0) ** 2 * 2.0)
    with pytest.raises(ValueError):
        verify_chi_margin_point(p, witness, 0.0)
    with pytest.raises(ValueError):
        verify_chi_margin_point(p, witness, 1.1)
