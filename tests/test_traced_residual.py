import numpy as np
import pytest

from openai_ns_reconstruction.traced_residual import (
    formal_force_from_full_spacetime_jets,
    traced_residual_from_full_spacetime_jets,
)


def endpoint_jets(degree: int, x: float, y: float, z: float) -> np.ndarray:
    if degree == 0:
        return np.array([1.0 + x, 2.0 + y, 3.0 + z])
    tensor = np.zeros((3,) + (4,) * degree)
    if degree == 1:
        tensor[:, 0] = np.array([4.0, -2.0, 1.0])
    return tensor


def test_traced_residual_uses_exact_extend_trace_branch() -> None:
    calls = {"past": 0, "jets": 0}

    def past(x: float, y: float, z: float, t: float) -> np.ndarray:
        calls["past"] += 1
        if t >= 1.0:
            raise AssertionError("past branch must not be sampled at/after endpoint")
        return np.array([t + x, y - t, z + 2.0 * t])

    def jets(degree: int, x: float, y: float, z: float) -> np.ndarray:
        calls["jets"] += 1
        return endpoint_jets(degree, x, y, z)

    traced = traced_residual_from_full_spacetime_jets(past, jets)

    np.testing.assert_allclose(traced(0.2, -0.3, 0.4, 0.75), [0.95, -1.05, 1.9])
    assert calls == {"past": 1, "jets": 0}

    np.testing.assert_allclose(traced(0.2, -0.3, 0.4, 1.0), [1.2, 1.7, 3.4])
    np.testing.assert_allclose(traced(0.2, -0.3, 0.4, 1.5), [1.2, 1.7, 3.4])
    assert calls == {"past": 1, "jets": 2}


def test_formal_force_glues_trace_to_normal_jet_borel_branch() -> None:
    def past(x: float, y: float, z: float, t: float) -> np.ndarray:
        return np.array([9.0 + t, 8.0 + x, 7.0 + z])

    force = formal_force_from_full_spacetime_jets(
        past,
        endpoint_jets,
        lambda degree: 1 << degree,
    )

    np.testing.assert_allclose(force(0.0, 0.0, 0.0, 0.5), [9.5, 8.0, 7.0])
    np.testing.assert_allclose(force(0.0, 0.0, 0.0, 1.0), [1.0, 2.0, 3.0])

    # At s=t-1=1/4, scale_0*s=1/4 and scale_1*s=1/2 are both on
    # the cutoff plateau, while scale_2*s=1 is already exactly zero.
    expected = np.array([1.0, 2.0, 3.0]) + 0.25 * np.array([4.0, -2.0, 1.0])
    np.testing.assert_allclose(force(0.0, 0.0, 0.0, 1.25), expected)

    # The degree-zero scale is at least one, hence the future branch is
    # structurally zero once t-endpoint >= 1.
    np.testing.assert_allclose(force(0.0, 0.0, 0.0, 2.0), np.zeros(3))


def test_traced_residual_fails_closed_on_bad_endpoint_data() -> None:
    def past(x: float, y: float, z: float, t: float) -> np.ndarray:
        return np.zeros(3)

    bad_shape = traced_residual_from_full_spacetime_jets(
        past, lambda degree, x, y, z: np.zeros((3, 4))
    )
    with pytest.raises(ValueError, match="order-0 full spacetime jet"):
        bad_shape(0.0, 0.0, 0.0, 1.0)

    nonfinite = traced_residual_from_full_spacetime_jets(
        past, lambda degree, x, y, z: np.array([np.nan, 0.0, 0.0])
    )
    with pytest.raises(ValueError, match="finite"):
        nonfinite(0.0, 0.0, 0.0, 1.0)

    with pytest.raises(ValueError, match="past_residual"):
        traced_residual_from_full_spacetime_jets(None, endpoint_jets)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="full_jets"):
        traced_residual_from_full_spacetime_jets(past, None)  # type: ignore[arg-type]
