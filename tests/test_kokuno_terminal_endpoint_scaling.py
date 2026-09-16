from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_terminal_endpoint_scaling import (
    EndpointScalingSignature,
    audit_fixed_period_identity,
    compose_signatures,
    fixed_period_viscosity_signature,
    parabolic_signature,
    viscosity_spatial_signature,
)


def test_viscosity_spatial_signature_exact_fixture() -> None:
    got = viscosity_spatial_signature(nu=Fraction(9, 16), sqrt_nu=Fraction(3, 4))
    assert got.velocity == Fraction(3, 4)
    assert got.pressure == Fraction(9, 16)
    assert got.force == Fraction(3, 4)
    assert got.x_argument == Fraction(4, 3)
    assert got.momentum_residual == Fraction(3, 4)
    assert got.divergence == 1
    assert got.energy_sq == Fraction(243, 1024)
    assert got.period == Fraction(3, 4)
    assert got.singular_time == 1


def test_parabolic_signature_exact_fixture() -> None:
    got = parabolic_signature(lam=Fraction(5, 3))
    assert got.velocity == Fraction(5, 3)
    assert got.pressure == Fraction(25, 9)
    assert got.force == Fraction(125, 27)
    assert got.x_argument == Fraction(5, 3)
    assert got.t_argument == Fraction(25, 9)
    assert got.momentum_residual == Fraction(125, 27)
    assert got.divergence == Fraction(25, 9)
    assert got.energy_sq == Fraction(3, 5)
    assert got.period == Fraction(3, 5)
    assert got.singular_time == Fraction(9, 25)


def test_fixed_period_identity_closes_exactly() -> None:
    composed, direct = audit_fixed_period_identity(nu=Fraction(9, 16), sqrt_nu=Fraction(3, 4))
    assert composed == direct
    assert direct == EndpointScalingSignature(
        velocity=Fraction(9, 16),
        pressure=Fraction(81, 256),
        force=Fraction(81, 256),
        x_argument=Fraction(1, 1),
        t_argument=Fraction(9, 16),
        momentum_residual=Fraction(81, 256),
        divergence=Fraction(9, 16),
        energy_sq=Fraction(81, 256),
        period=Fraction(1, 1),
        singular_time=Fraction(16, 9),
    )


def test_period_and_energy_factors_cancel_through_composition() -> None:
    spatial = viscosity_spatial_signature(nu=Fraction(25, 49), sqrt_nu=Fraction(5, 7))
    para = parabolic_signature(lam=Fraction(5, 7))
    composed = compose_signatures(spatial, para)
    assert composed.period == 1
    assert composed.energy_sq == Fraction(625, 2401)
    assert composed.singular_time == Fraction(49, 25)


def test_mutated_square_root_fails_closed() -> None:
    with pytest.raises(ValueError, match="square exactly"):
        viscosity_spatial_signature(
            nu=Fraction(9, 16),
            sqrt_nu=Fraction(3, 4) + Fraction(1, 2**40),
        )


def test_approximate_and_nonpositive_inputs_fail_closed() -> None:
    with pytest.raises(TypeError):
        fixed_period_viscosity_signature(nu=0.5)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        parabolic_signature(lam=1.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        fixed_period_viscosity_signature(nu=Fraction(0, 1))
    with pytest.raises(ValueError):
        parabolic_signature(lam=Fraction(-1, 2))


def test_composition_rejects_non_signature() -> None:
    sig = fixed_period_viscosity_signature(nu=Fraction(1, 1))
    with pytest.raises(TypeError):
        compose_signatures(sig, object())  # type: ignore[arg-type]
