from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_radial_tensor_force import (
    FULL_RECONSTRUCTION,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    cartesian_tensor_divergence,
    cylindrical_force_components,
    cylindrical_to_cartesian,
    radial_inverse_rhs,
)

R = F(3, 2)
S1 = F(5, 7)
DR1 = F(-2, 9)
DZ1 = F(11, 13)
S2 = F(-4, 5)
DR2 = F(7, 8)


def test_public_r34_components_exact():
    assert cylindrical_force_components(R, S1, DR1, DZ1, S2, DR2) == (
        F(11, 13),
        F(-23, 120),
        F(16, 63),
    )


def test_radial_inverse_binding_exact():
    f2 = F(59, 120)
    b2m2 = F(3, 10)
    f1 = F(-53, 126)
    b1m1 = F(-1, 6)
    force = cylindrical_force_components(R, S1, DR1, DZ1, S2, DR2)
    assert radial_inverse_rhs(f2, b2m2) == force[1]
    assert radial_inverse_rhs(f1, b1m1) == force[2]


def test_direct_cartesian_product_rule_matches_cylindrical_formula():
    x, y = F(9, 10), F(6, 5)
    direct = cartesian_tensor_divergence(x, y, R, S1, DR1, DZ1, S2, DR2)
    expected = cylindrical_to_cartesian(
        cylindrical_force_components(R, S1, DR1, DZ1, S2, DR2), x, y, R
    )
    assert direct == expected


def test_omitting_one_cylindrical_connection_term_is_detected():
    corrected_theta = cylindrical_force_components(R, S1, DR1, DZ1, S2, DR2)[1]
    wrong_theta = DR2 + S2 / R
    assert corrected_theta - wrong_theta == S2 / R == F(-8, 15)


def test_dropping_radial_force_is_detected():
    force = cylindrical_force_components(R, S1, DR1, DZ1, S2, DR2)
    assert force[0] == F(11, 13)
    assert force[0] != 0


def test_fail_closed_axis_and_type_boundaries():
    with pytest.raises(ValueError):
        cylindrical_force_components(0, S1, DR1, DZ1, S2, DR2)
    with pytest.raises(TypeError):
        cylindrical_force_components(1.0, S1, DR1, DZ1, S2, DR2)


def test_truth_flags_stay_false():
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
