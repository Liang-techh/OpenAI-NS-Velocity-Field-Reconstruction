from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section10_paper_time_switch_endpoint_jet import (
    Section10PaperTimeSwitchEndpointJetCertificate,
)


def _scalar_jet(order: int) -> dict[tuple[int, int, int, int], float]:
    return {
        alpha: float(1 + 2 * alpha[0] + 3 * alpha[1] + 5 * alpha[2] + 7 * alpha[3])
        for alpha in required_spacetime_multiindices(order)
    }


def _vector_jet(order: int) -> dict[tuple[int, int, int, int], np.ndarray]:
    return {
        alpha: np.array(
            [
                1.0 + alpha[0] + 2 * alpha[1],
                2.0 + 3 * alpha[2] + alpha[3],
                -4.0 + alpha[0] + alpha[1] + alpha[2] + alpha[3],
            ]
        )
        for alpha in required_spacetime_multiindices(order)
    }


def test_pinned_endpoint_switch_jet_is_exact_identity_multiplier() -> None:
    certificate = Section10PaperTimeSwitchEndpointJetCertificate.pinned(4)

    assert certificate.endpoint == Fraction(1, 1)
    assert certificate.switch_time_jet == (
        Fraction(1, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    )
    assert certificate.endpoint_multiplier_jet_is_identity is True


def test_scalar_endpoint_jet_uses_leibniz_and_preserves_mixed_derivatives() -> None:
    certificate = Section10PaperTimeSwitchEndpointJetCertificate.pinned(3)
    source = _scalar_jet(3)
    switched = certificate.apply_scalar_spacetime_jet(source)

    assert switched == source

    # Independent hand expansion for D_t^2 D_x(sigma F) at t=1:
    # sigma F_ttx + 2 sigma_t F_tx + sigma_tt F_x.
    alpha = (2, 1, 0, 0)
    manual = (
        source[(2, 1, 0, 0)]
        + 2.0 * 0.0 * source[(1, 1, 0, 0)]
        + 0.0 * source[(0, 1, 0, 0)]
    )
    assert switched[alpha] == manual


def test_vector_endpoint_jet_is_transparent_componentwise() -> None:
    certificate = Section10PaperTimeSwitchEndpointJetCertificate.pinned(2)
    source = _vector_jet(2)
    switched = certificate.apply_vector_spacetime_jet(source)

    for alpha in required_spacetime_multiindices(2):
        np.testing.assert_array_equal(switched[alpha], source[alpha])
        assert switched[alpha].flags.writeable is False


def test_endpoint_adapter_fails_closed_on_incomplete_or_nonfinite_jet() -> None:
    certificate = Section10PaperTimeSwitchEndpointJetCertificate.pinned(2)
    source = _scalar_jet(2)
    source.pop((0, 0, 0, 0))

    with pytest.raises(ValueError, match="cover exactly"):
        certificate.apply_scalar_spacetime_jet(source)

    vector = _vector_jet(2)
    vector[(1, 0, 0, 0)] = np.array([1.0, np.nan, 3.0])
    with pytest.raises(ValueError, match="finite three-vector"):
        certificate.apply_vector_spacetime_jet(vector)


def test_certificate_does_not_claim_missing_t1_extension_or_final_velocity() -> None:
    certificate = Section10PaperTimeSwitchEndpointJetCertificate.pinned(5)

    assert certificate.section9_endpoint_jet_supplied_by_actual_construction is False
    assert certificate.section9_field_smooth_extension_through_t1_constructed is False
    assert certificate.endpoint_residual_closure_verified is False
    assert certificate.paper_exact_velocity_available is False

    with pytest.raises(ValueError, match="nonnegative integer"):
        Section10PaperTimeSwitchEndpointJetCertificate.pinned(-1)
    with pytest.raises(ValueError, match="exceeds certificate order"):
        certificate.switch_time_derivative(6)
