from dataclasses import dataclass
from fractions import Fraction as Q

import pytest

from openai_ns_reconstruction.kokuno_profile_moment_scaling import (
    FULL_RECONSTRUCTION,
    IMPORTED_PROFILE_EXISTENCE_PROVED_HERE,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    FiveMomentState,
    FixedRectangleScale,
    bind_profile_moments,
    from_fixed_rectangle,
    to_fixed_rectangle,
)


def test_exact_five_moment_scaling_and_inverse() -> None:
    hat = FiveMomentState(M=Q(2, 3), I=Q(-5, 7), J=Q(11, 13), S=Q(17, 19), C_p=Q(-23, 29))
    scale = FixedRectangleScale(Q(3))  # rho=9, rho^(3/2)=27 exactly

    physical = from_fixed_rectangle(hat, scale=scale)
    assert physical == FiveMomentState(
        M=Q(6), I=Q(-135, 7), J=Q(297, 13), S=Q(153, 19), C_p=Q(-23, 29)
    )
    assert to_fixed_rectangle(physical, scale=scale) == hat


def test_all_five_components_participate_in_exact_match() -> None:
    baseline = FiveMomentState(M=1, I=2, J=3, S=4, C_p=5)
    assert baseline.residual_against(baseline).is_zero

    names = ("M", "I", "J", "S", "C_p")
    for name in names:
        values = dict(M=1, I=2, J=3, S=4, C_p=5)
        values[name] = Q(values[name]) + Q(1, 2**40)
        perturbed = FiveMomentState(**values)
        residual = perturbed.residual_against(baseline)
        assert not residual.is_zero
        assert getattr(residual, name) == Q(1, 2**40)


def test_exact_api_rejects_float_scales_and_values() -> None:
    with pytest.raises(TypeError):
        FixedRectangleScale(3.0)
    with pytest.raises(TypeError):
        FiveMomentState(M=1.0, I=2, J=3, S=4, C_p=5)
    with pytest.raises(ValueError):
        FixedRectangleScale(0)


@dataclass
class DummyLeadingProfile:
    name: str = "typed-profile-fixture"
    paper_exact: bool = False
    provenance: str | None = "independent test fixture"


def test_binding_accepts_target_leading_profile_shape_but_requires_provenance() -> None:
    moments = FiveMomentState(M=1, I=2, J=3, S=4, C_p=5)
    bound = bind_profile_moments(
        DummyLeadingProfile(), X_h=Q(7, 5), moments=moments, data_provenance="fixture-five-moments"
    )
    assert bound.profile_name == "typed-profile-fixture"
    assert bound.X_h == Q(7, 5)
    assert bound.moments == moments

    with pytest.raises(ValueError):
        bind_profile_moments(
            DummyLeadingProfile(provenance=None),
            X_h=1,
            moments=moments,
            data_provenance="fixture-five-moments",
        )
    with pytest.raises(ValueError):
        bind_profile_moments(
            DummyLeadingProfile(), X_h=1, moments=moments, data_provenance=""
        )


def test_truth_boundary_stays_fail_closed() -> None:
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
    assert IMPORTED_PROFILE_EXISTENCE_PROVED_HERE is False
