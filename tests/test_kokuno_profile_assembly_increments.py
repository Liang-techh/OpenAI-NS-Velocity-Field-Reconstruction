from decimal import Decimal
from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_profile_assembly_increments import (
    FiveMomentDensity,
    ProfileAssemblyInput,
    all_zero,
    direct_density_difference,
    pa14_increment_formula,
    pa14_residual,
)


def fixture() -> ProfileAssemblyInput:
    return ProfileAssemblyInput.exact(
        X=F(9, 32),
        radial_factor=F(3, 4),
        U=F(5, 7),
        E=F(-4, 9),
        u=F(2, 11),
        e=F(3, 13),
    )


def test_exact_pa14_formula_matches_direct_density_difference() -> None:
    data = fixture()
    expected = FiveMomentDensity(
        M=F(2, 11),
        I=F(9, 52),
        J=F(1135, 12012),
        S=F(316691, 858858),
        C_p=F(-1232, 4563),
    )
    assert pa14_increment_formula(data) == expected
    assert direct_density_difference(data) == expected
    assert all_zero(pa14_residual(data))


def test_cross_term_in_J_is_not_optional() -> None:
    data = fixture()
    exact = pa14_increment_formula(data)
    H = data.radial_factor * data.E
    without_cross = H * data.u + data.U * data.radial_factor * data.e
    assert exact.J - without_cross == F(9, 286)
    assert exact.J != without_cross


def test_quadratic_terms_in_S_are_not_optional() -> None:
    data = fixture()
    exact = pa14_increment_formula(data)
    linearized = 2 * data.U * data.u - data.E * data.e
    assert exact.S - linearized == F(263, 40898)
    assert exact.S != linearized


def test_quadratic_term_in_Cp_is_not_optional() -> None:
    data = fixture()
    exact = pa14_increment_formula(data)
    linearized = data.E * data.e / data.X
    assert exact.C_p - linearized == F(16, 169)
    assert exact.C_p != linearized


def test_radial_factor_mutation_fails_closed() -> None:
    with pytest.raises(ValueError, match=r"radial_factor\*\*2"):
        ProfileAssemblyInput.exact(
            X=F(9, 32),
            radial_factor=F(3, 4) + F(1, 2**40),
            U=F(5, 7),
            E=F(-4, 9),
            u=F(2, 11),
            e=F(3, 13),
        )


def test_approximate_theorem_inputs_are_rejected() -> None:
    with pytest.raises(TypeError):
        ProfileAssemblyInput.exact(
            X=0.28125, radial_factor=F(3, 4), U=0, E=0, u=0, e=0
        )
    with pytest.raises(TypeError):
        ProfileAssemblyInput.exact(
            X=F(9, 32), radial_factor=F(3, 4), U=Decimal("0.5"), E=0, u=0, e=0
        )


def test_domain_boundary_rejects_nonpositive_X_and_boolean_inputs() -> None:
    with pytest.raises(ValueError):
        ProfileAssemblyInput.exact(X=0, radial_factor=1, U=0, E=0, u=0, e=0)
    with pytest.raises(TypeError):
        ProfileAssemblyInput.exact(
            X=F(9, 32), radial_factor=F(3, 4), U=True, E=0, u=0, e=0
        )
