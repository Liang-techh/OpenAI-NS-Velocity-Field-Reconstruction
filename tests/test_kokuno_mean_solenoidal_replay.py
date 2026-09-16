from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_mean_solenoidal_replay import MeanSolenoidalJet


def _fixture() -> MeanSolenoidalJet:
    return MeanSolenoidalJet(
        radius=Fraction(3, 2),
        epsilon=Fraction(1, 100),
        psi=Fraction(7, 5),
        d_r_psi=Fraction(-11, 7),
        d_z_psi=Fraction(13, 17),
        d_rz_psi=Fraction(-19, 23),
    )


def test_mc14_solenoidal_identity_is_exact():
    jet = _fixture()
    assert jet.delta_beta == Fraction(-13, 1700)
    assert jet.delta_gamma == Fraction(-67, 105)
    assert jet.d_r_delta_beta == Fraction(19, 2300)
    assert jet.d_z_delta_gamma == Fraction(-371, 1173)
    assert jet.cylindrical_divergence == 0
    assert jet.identity_holds_exactly is True


def test_required_minus_sign_fails_closed_when_flipped():
    jet = _fixture()
    assert jet.divergence_with_beta_sign(-1) == 0
    assert jet.divergence_with_beta_sign(1) == Fraction(-371, 58650)
    assert jet.divergence_with_beta_sign(1) != 0


def test_approximate_inputs_are_rejected():
    with pytest.raises(TypeError):
        MeanSolenoidalJet(
            radius=1.0,
            epsilon=Fraction(1, 100),
            psi=0,
            d_r_psi=0,
            d_z_psi=0,
            d_rz_psi=0,
        )
    with pytest.raises(TypeError):
        MeanSolenoidalJet(
            radius=1,
            epsilon=Decimal("0.01"),
            psi=0,
            d_r_psi=0,
            d_z_psi=0,
            d_rz_psi=0,
        )


def test_domain_and_sign_controls_are_fail_closed():
    with pytest.raises(ValueError):
        MeanSolenoidalJet(
            radius=0,
            epsilon=1,
            psi=0,
            d_r_psi=0,
            d_z_psi=0,
            d_rz_psi=0,
        )
    with pytest.raises(ValueError):
        MeanSolenoidalJet(
            radius=1,
            epsilon=0,
            psi=0,
            d_r_psi=0,
            d_z_psi=0,
            d_rz_psi=0,
        )
    with pytest.raises(ValueError):
        _fixture().divergence_with_beta_sign(0)


def test_truth_boundary_remains_fail_closed():
    jet = _fixture()
    assert jet.status == "formula-level-replay"
    assert jet.actual_compact_correction_field_materialized is False
    assert jet.actual_wave_defect_consumed is False
    assert jet.all_stage_convergence_certified is False
    assert jet.paper_exact_velocity_available is False
    assert jet.full_reconstruction is False
