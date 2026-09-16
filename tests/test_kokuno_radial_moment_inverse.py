from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_radial_moment_inverse import (
    normalize_profile,
    replay_radial_moment_inverse,
    weighted_moment,
)


def fixture(e: int):
    lower = F(1, 2)
    upper = F(5, 2)
    radius = F(7, 5)
    # Independent exact polynomial fixtures; no paper stress/profile existence is asserted.
    source = (F(7, 9), F(-5, 11), F(2, 7), F(1, 13))
    raw_profile = (F(3, 5), F(1, 7), F(2, 9))
    profile = normalize_profile(raw_profile, e=e, lower=lower, upper=upper)
    return lower, upper, radius, source, profile


@pytest.mark.parametrize("e", [1, 2])
def test_exact_moment_complement_and_radial_inverse(e):
    lower, upper, radius, source, profile = fixture(e)
    result = replay_radial_moment_inverse(
        e=e,
        source=source,
        normalized_profile=profile,
        lower=lower,
        upper=upper,
        radius=radius,
    )
    assert result.profile_moment == 1
    assert result.projected_moment == 0
    assert result.upper_weighted_primitive == 0
    assert result.inverse_residual == 0
    assert result.verified


def test_nonzero_source_moment_is_retained_not_discarded():
    lower, upper, radius, source, profile = fixture(2)
    result = replay_radial_moment_inverse(
        e=2,
        source=source,
        normalized_profile=profile,
        lower=lower,
        upper=upper,
        radius=radius,
    )
    assert result.source_moment != 0
    assert result.expected_differential_value == -result.projected_source_value


def test_profile_normalization_drift_fails_closed():
    lower, upper, radius, source, profile = fixture(1)
    drifted = list(profile)
    drifted[0] += F(1, 2**40)
    assert weighted_moment(drifted, e=1, lower=lower, upper=upper) != 1
    with pytest.raises(ValueError, match="weighted moment 1"):
        replay_radial_moment_inverse(
            e=1,
            source=source,
            normalized_profile=drifted,
            lower=lower,
            upper=upper,
            radius=radius,
        )


def test_approximate_inputs_rejected():
    lower, upper, radius, source, profile = fixture(1)
    with pytest.raises(TypeError):
        replay_radial_moment_inverse(
            e=1,
            source=(0.1,),
            normalized_profile=profile,
            lower=lower,
            upper=upper,
            radius=radius,
        )


def test_domain_and_weight_restrictions_fail_closed():
    lower, upper, radius, source, profile = fixture(1)
    with pytest.raises(ValueError, match="e must be 1 or 2"):
        replay_radial_moment_inverse(
            e=3,
            source=source,
            normalized_profile=profile,
            lower=lower,
            upper=upper,
            radius=radius,
        )
    with pytest.raises(ValueError, match="radius must lie"):
        replay_radial_moment_inverse(
            e=1,
            source=source,
            normalized_profile=profile,
            lower=lower,
            upper=upper,
            radius=F(3),
        )
