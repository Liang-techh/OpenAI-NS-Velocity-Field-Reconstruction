from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_base_heat_coordinate_bridge import (
    exact_fraction_power,
    verify_base_heat_coordinate_bridge,
)


def fixture():
    h = Fraction(1, 200)
    q = Fraction(3, 2) ** 200
    X = Fraction(5, 4) ** 200
    eta = Fraction(1, 3)
    d = 1 - eta * eta
    tau = q * d
    s = q * X
    return dict(
        q=q,
        X=X,
        eta=eta,
        d=d,
        tau=tau,
        s=s,
        h=h,
        c_inf=Fraction(17, 19),
        H_value=Fraction(11, 13),
    )


def test_exact_coordinate_and_heat_value_identity():
    data = fixture()
    result = verify_base_heat_coordinate_bridge(**data)
    assert result.A == Fraction(101, 200)
    assert result.z_profile == Fraction(16, 9) * Fraction(4, 5) ** 200
    assert result.z_physical == result.z_profile
    expected = Fraction(17, 19) * Fraction(11, 13) * Fraction(8, 15) ** 101
    assert result.profile_value == expected
    assert result.physical_value == expected
    assert result.z_residual == 0
    assert result.value_residual == 0


def test_exact_fraction_power_requires_certified_root():
    assert exact_fraction_power(Fraction(3, 2) ** 200, Fraction(-101, 200)) == Fraction(2, 3) ** 101
    with pytest.raises(ValueError):
        exact_fraction_power(Fraction(2), Fraction(-101, 200))


def test_wrong_tau_fails_closed_at_exact_2_to_minus_40_mutation():
    data = fixture()
    data["tau"] += Fraction(1, 2**40)
    with pytest.raises(ValueError, match=r"tau must equal q\*d exactly"):
        verify_base_heat_coordinate_bridge(**data)


def test_wrong_s_convention_fails_closed():
    data = fixture()
    data["s"] *= 2
    with pytest.raises(ValueError, match=r"s=q\*X exactly"):
        verify_base_heat_coordinate_bridge(**data)


def test_wrong_d_fails_closed():
    data = fixture()
    data["d"] += Fraction(1, 2**40)
    with pytest.raises(ValueError, match=r"1-eta\^2"):
        verify_base_heat_coordinate_bridge(**data)


def test_float_theorem_input_is_rejected():
    data = fixture()
    data["h"] = 0.005
    with pytest.raises(TypeError):
        verify_base_heat_coordinate_bridge(**data)


def test_source_h_boundary_is_rejected():
    data = fixture()
    data["h"] = Fraction(1, 100)
    with pytest.raises(ValueError, match="0 < h < 1/100"):
        verify_base_heat_coordinate_bridge(**data)
