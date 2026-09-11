import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_support import (
    formal_force_support_point_certificate,
)
from openai_ns_reconstruction.spatial_localization import support_cylinder_contains


OUTSIDE = (0.30, 0.0, 0.0)


def zero_past(x, y, z, t):
    assert not support_cylinder_contains(x, y, z)
    return np.zeros(3)


def supported_full_jets(degree, x, y, z):
    tensor = np.zeros((3,) + (4,) * degree)
    if support_cylinder_contains(x, y, z):
        tensor[(0,) + (0,) * degree] = degree + 1.0
    return tensor


def unit_local_scale(_degree):
    return 1


def test_past_branch_does_not_touch_endpoint_data() -> None:
    def unavailable_jets(*_args):
        raise AssertionError("past branch must not request endpoint jets")

    certificate = formal_force_support_point_certificate(
        zero_past,
        unavailable_jets,
        unit_local_scale,
        *OUTSIDE,
        0.9,
    )
    assert certificate.branch == "past"
    assert certificate.checked_degrees == ()
    assert certificate.certified


def test_endpoint_checks_only_degree_zero_trace() -> None:
    calls = []

    def jets(degree, x, y, z):
        calls.append(degree)
        if degree != 0:
            raise AssertionError("t=T needs only the degree-zero coefficient")
        return np.zeros(3)

    certificate = formal_force_support_point_certificate(
        zero_past,
        jets,
        unit_local_scale,
        *OUTSIDE,
        1.0,
    )
    assert certificate.branch == "endpoint"
    assert certificate.checked_degrees == (0,)
    assert certificate.certified
    assert calls and set(calls) == {0}


def test_future_certificate_uses_only_locally_finite_active_coefficients() -> None:
    calls = []

    def jets(degree, x, y, z):
        calls.append(degree)
        if degree > 2:
            raise AssertionError("inactive Taylor coefficients must not be sampled")
        return np.zeros((3,) + (4,) * degree)

    # With local scale 1 the doubling envelope is 1,2,4,8,... .  At s=0.2,
    # degrees 0,1,2 can contribute and degree 3 is already outside cutoff support.
    certificate = formal_force_support_point_certificate(
        zero_past,
        jets,
        unit_local_scale,
        *OUTSIDE,
        1.2,
    )
    assert certificate.branch == "future"
    assert certificate.checked_degrees == (0, 1, 2)
    assert certificate.certified
    assert calls and max(calls) == 2


def test_nonzero_active_coefficient_fails_closed() -> None:
    def jets(degree, x, y, z):
        tensor = np.zeros((3,) + (4,) * degree)
        if degree == 1:
            tensor[(0, 0)] = 3.0
        return tensor

    certificate = formal_force_support_point_certificate(
        zero_past,
        jets,
        unit_local_scale,
        *OUTSIDE,
        1.2,
    )
    assert certificate.checked_degrees == (0, 1, 2)
    assert not certificate.relevant_inputs_zero
    assert not certificate.formal_value_zero
    assert not certificate.certified


def test_nonzero_past_residual_fails_closed() -> None:
    def nonzero_past(_x, _y, _z, _t):
        return np.array([1.0, 0.0, 0.0])

    certificate = formal_force_support_point_certificate(
        nonzero_past,
        supported_full_jets,
        unit_local_scale,
        *OUTSIDE,
        0.9,
    )
    assert not certificate.certified


def test_certificate_refuses_points_inside_closed_support_cylinder() -> None:
    with pytest.raises(ValueError, match="outside the closed Section 10 cylinder"):
        formal_force_support_point_certificate(
            zero_past,
            supported_full_jets,
            unit_local_scale,
            0.0,
            0.0,
            0.0,
            1.2,
        )
