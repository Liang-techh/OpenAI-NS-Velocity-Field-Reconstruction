import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_limit_majorant import EndpointPowerLawMajorant
from openai_ns_reconstruction.section10_endpoint_borel_admission import (
    admit_section10_endpoint_candidate_for_borel_prefix,
)
from openai_ns_reconstruction.section10_endpoint_ladder import Section10EndpointMajorantLadder


NORMAL_JETS = (
    np.array([0.8, -0.2, 0.35]),
    np.array([-0.15, 0.4, 0.1]),
    np.array([0.3, 0.05, -0.25]),
)


def _endpoint_full_jets(degree: int, x: float, y: float, z: float) -> np.ndarray:
    shape = (3,) + (4,) * degree
    out = np.zeros(shape)
    normal = NORMAL_JETS[degree].copy()
    if degree == 0:
        normal += np.array([0.02 * x, -0.03 * y, 0.01 * z])
    out[(slice(None),) + (0,) * degree] = normal
    # Put a non-normal component in higher jets so the left trace gate really
    # checks the full dense tensor rather than only the future normal jet.
    if degree >= 1:
        out[(0,) + (1,) * degree] = 0.07 * degree
    return out


def _drift_tensor(degree: int) -> np.ndarray:
    out = np.zeros((3,) + (4,) * degree)
    out.reshape(-1)[-1] = 0.1
    return out


def _past_full_jets(
    degree: int, x: float, y: float, z: float, t: float
) -> np.ndarray:
    # Analytic fixture with true endpoint _endpoint_full_jets and constant time
    # derivative norm 0.1.  The independently supplied coefficient 0.2 below
    # dominates that derivative, so its alpha=0 tail radius is 0.2*(1-t).
    return _endpoint_full_jets(degree, x, y, z) + (1.0 - t) * _drift_tensor(degree)


def _ladder() -> Section10EndpointMajorantLadder:
    return Section10EndpointMajorantLadder(
        EndpointPowerLawMajorant(
            derivative_degree=degree,
            spatial_window=3,
            coefficient=0.2,
            singularity_exponent=0.0,
            valid_from=0.75,
            endpoint=1.0,
        )
        for degree in range(3)
    )


def test_true_analytic_endpoint_is_admitted_then_frozen_for_right_prefix() -> None:
    point = (0.11, -0.08, 0.04)
    calls = {0: 0, 1: 0, 2: 0}

    def counted_candidate(degree: int, x: float, y: float, z: float) -> np.ndarray:
        calls[degree] += 1
        return _endpoint_full_jets(degree, x, y, z)

    certificate = admit_section10_endpoint_candidate_for_borel_prefix(
        _ladder(),
        _past_full_jets,
        counted_candidate,
        lambda degree: 1 << degree,
        *point,
        0.9,
    )

    assert certificate.formal_candidate_ready
    assert certificate.right_prefix.formal_prefix_ready
    assert calls == {0: 1, 1: 1, 2: 1}

    # Independent closed-form check for this fixture: the center is exactly
    # 0.01 away from the true endpoint in one dense coordinate, while the
    # alpha=0 majorant gives radius 0.02 at t=0.9.
    for row in certificate.rows:
        assert row.distance_to_center == pytest.approx(0.01, abs=2e-15)
        assert row.enclosure_radius == pytest.approx(0.02, abs=2e-15)
        assert row.admitted

    # The right prefix consumes exactly the frozen endpoint candidate and its
    # normal time contractions agree with the known analytic endpoint tensors.
    for degree, row in enumerate(certificate.right_prefix.rows):
        expected = _endpoint_full_jets(degree, *point)[
            (slice(None),) + (0,) * degree
        ]
        assert row.prescribed_normal_jet == pytest.approx(expected, rel=0.0, abs=0.0)
        assert row.prefix_endpoint_time_jet == pytest.approx(expected, rel=0.0, abs=0.0)


def test_candidate_outside_one_trace_ball_fails_before_borel_scale_is_used() -> None:
    point = (0.03, 0.02, -0.01)

    def bad_candidate(degree: int, x: float, y: float, z: float) -> np.ndarray:
        out = _endpoint_full_jets(degree, x, y, z)
        if degree == 1:
            out = out.copy()
            out.reshape(-1)[0] += 0.05
        return out

    def forbidden_scale(degree: int) -> int:
        raise AssertionError("Borel scale must not be queried after left-trace rejection")

    with pytest.raises(ValueError, match="degree-1 endpoint candidate lies outside"):
        admit_section10_endpoint_candidate_for_borel_prefix(
            _ladder(),
            _past_full_jets,
            bad_candidate,
            forbidden_scale,
            *point,
            0.9,
        )


def test_malformed_candidate_dense_jet_fails_closed() -> None:
    def wrong_shape(degree: int, x: float, y: float, z: float) -> np.ndarray:
        if degree == 2:
            return np.zeros((3, 4))
        return _endpoint_full_jets(degree, x, y, z)

    with pytest.raises(ValueError, match="order-2 full spacetime jet"):
        admit_section10_endpoint_candidate_for_borel_prefix(
            _ladder(),
            _past_full_jets,
            wrong_shape,
            lambda degree: 1 << degree,
            0.0,
            0.0,
            0.0,
            0.9,
        )


def test_admission_keeps_actual_limit_and_infinite_smoothness_boundaries_closed() -> None:
    certificate = admit_section10_endpoint_candidate_for_borel_prefix(
        _ladder(),
        _past_full_jets,
        _endpoint_full_jets,
        lambda degree: 2 * (1 << degree),
        0.0,
        0.0,
        0.0,
        0.95,
    )

    assert certificate.endpoint_trace_candidate_consistency_verified
    assert certificate.finite_prefix_right_jet_algebra_verified
    assert not certificate.source_majorants_derived_from_actual_residual_verified
    assert not certificate.actual_section9_residual_limits_verified
    assert not certificate.endpoint_limit_uniqueness_verified
    assert not certificate.analytic_template_bounds_verified
    assert not certificate.infinite_borel_right_jets_verified
    assert not certificate.all_order_borel_smoothness_verified
    assert not certificate.smooth_compact_forcing_verified
    assert not certificate.paper_exact_velocity_available
