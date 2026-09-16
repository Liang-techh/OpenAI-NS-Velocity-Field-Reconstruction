from fractions import Fraction
import math
import pytest

from openai_ns_reconstruction.kokuno_profile_radial_primitive import (
    FULL_RECONSTRUCTION,
    IMPORTED_PROFILE_EXISTENCE_PROVED_HERE,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    RadialPrimitiveJet,
    compare_target_leading_profile,
)


def _fixture(**changes):
    X = Fraction(1, 2)
    eta = Fraction(1, 2)
    h = Fraction(1, 16)
    U = 1 + X / 2 + eta / 4
    U_X = Fraction(1, 2)
    U_eta = Fraction(1, 4)
    M = X * (1 + X / 4 + eta / 4)
    M_eta = X * U_eta
    data = dict(X=X, eta=eta, h=h, U=U, U_X=U_X, U_eta=U_eta,
                M=M, M_eta=M_eta, M_X=U, M_eta_X=U_eta)
    data.update(changes)
    return RadialPrimitiveJet(**data)


def test_exact_radial_primitive_closes_divergence():
    jet = _fixture()
    assert jet.primitive_compatibility_residuals == (0, 0)
    assert jet.V0() == Fraction(41, 124)
    assert jet.dV0_dX_from_primitive() == Fraction(107, 124)
    assert jet.divergence_rhs() == Fraction(107, 124)
    assert jet.divergence_residual() == 0


def test_primitive_derivative_defect_is_not_tolerated():
    good = _fixture()
    bad = _fixture(M_X=good.U + Fraction(1, 2**40))
    assert bad.primitive_compatibility_residuals[0] == Fraction(1, 2**40)
    assert bad.divergence_residual() != 0
    assert bad.divergence_residual() == -Fraction(7, 17_042_430_230_528)


def test_exact_inputs_fail_closed():
    good = _fixture()
    with pytest.raises(TypeError):
        _fixture(U=float(good.U))
    with pytest.raises(ValueError):
        _fixture(h=Fraction(1, 2))
    with pytest.raises(ValueError):
        _fixture(eta=Fraction(5, 4))


def test_target_leading_profile_bridge_matches_existing_typed_interface():
    from openai_ns_reconstruction.profiles import LeadingProfile

    def U(X, eta):
        return 1.0 + X / 2.0 + eta / 4.0

    profile = LeadingProfile(
        E=lambda X, eta: 0.0,
        U=U,
        dU_deta=lambda X, eta: 0.25,
        average_U=lambda X, eta: 1.0 + X / 4.0 + eta / 4.0,
        average_dU_deta=lambda X, eta: 0.25,
        name="kokuno-radial-primitive-regression-fixture",
        provenance="analytic linear-U fixture for Kokuno NS-profile radial primitive bridge",
    )
    result = compare_target_leading_profile(profile, _fixture())
    assert abs(result.V0_error) <= 2 * math.ulp(float(result.expected_V0))
    assert abs(result.flux_error) <= 2 * math.ulp(float(result.expected_flux_factor))


def test_truth_boundary_stays_fail_closed():
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
    assert IMPORTED_PROFILE_EXISTENCE_PROVED_HERE is False
