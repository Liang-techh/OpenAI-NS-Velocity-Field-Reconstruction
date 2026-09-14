import pytest

from openai_ns_reconstruction.profiles import LeadingProfile


def _profile(*, U, average_U=None):
    return LeadingProfile(
        E=lambda X, eta: 0.0,
        U=U,
        dU_deta=lambda X, eta: 0.0,
        average_U=average_U,
        paper_exact=False,
    )


def test_exact_average_callback_is_authoritative_at_axis():
    def forbidden_u(X: float, eta: float) -> float:
        raise AssertionError("exact average callback must not sample U")

    profile = _profile(
        U=forbidden_u,
        average_U=lambda X, eta: 7.0 + X - eta,
    )

    assert profile.radial_average_U(0.0, 0.25, n=2) == pytest.approx(6.75)


def test_axis_average_falls_back_to_profile_value_without_exact_callback():
    seen = []

    def U(X: float, eta: float) -> float:
        seen.append((X, eta))
        return 3.5 + eta

    profile = _profile(U=U)

    assert profile.radial_average_U(0.0, -0.25, n=2) == pytest.approx(3.25)
    assert seen == [(0.0, -0.25)]
