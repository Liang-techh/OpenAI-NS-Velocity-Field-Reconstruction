import numpy as np
import pytest

from openai_ns_reconstruction.cutoffs import standard_cutoff
from openai_ns_reconstruction.endpoint_borel import (
    BorelRightExtension,
    DoublingEnvelope,
    borel_monomial_coefficient,
)


def test_doubling_envelope_matches_pinned_recurrence() -> None:
    local = (0, 3, 2, 20, 1)
    envelope = DoublingEnvelope(lambda j: local[j] if j < len(local) else 0)
    assert [envelope(j) for j in range(6)] == [1, 3, 6, 20, 40, 80]


def test_right_extension_is_the_locally_finite_borel_sum() -> None:
    def jet(j: int, x: float, y: float, z: float) -> np.ndarray:
        return np.array([j + 1.0 + x, (-1.0) ** j + y, z + 0.25 * j])

    extension = BorelRightExtension(jet, lambda j: 1, endpoint=1.0)
    point = (0.2, -0.1, 0.3, 1.1)
    s = point[3] - 1.0

    # Independent finite sum: envelope is 1,2,4,8,16,... and the first term
    # with scale*s >= 1, together with all later terms, has zero cutoff.
    expected = np.zeros(3)
    scale = 1
    for j in range(20):
        if scale * s >= 1.0:
            break
        expected += (
            standard_cutoff(scale * s)
            * borel_monomial_coefficient(s, j)
            * jet(j, *point[:3])
        )
        scale *= 2

    assert np.allclose(extension(*point), expected, rtol=1e-14, atol=1e-14)


def test_right_boundary_jets_match_a_manufactured_finite_jet_family() -> None:
    a0 = np.array([1.2, -0.7, 0.4])
    a1 = np.array([-0.3, 0.8, 1.1])
    a2 = np.array([0.6, -1.4, 0.2])

    def jet(j: int, x: float, y: float, z: float) -> np.ndarray:
        if j == 0:
            return a0
        if j == 1:
            return a1
        if j == 2:
            return a2
        return np.zeros(3)

    extension = BorelRightExtension(jet, lambda j: 1, endpoint=1.0)
    h = 1e-4
    f0 = extension(0.0, 0.0, 0.0, 1.0)
    f1 = extension(0.0, 0.0, 0.0, 1.0 + h)
    f2 = extension(0.0, 0.0, 0.0, 1.0 + 2.0 * h)

    # For these three nonzero degrees, all cutoffs are exactly on their
    # plateau at 0,h,2h, so the one-sided formulas independently recover the
    # prescribed Taylor jets (up to floating arithmetic).
    first = (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
    second = (f0 - 2.0 * f1 + f2) / (h * h)
    assert np.allclose(f0, a0, rtol=0.0, atol=1e-14)
    assert np.allclose(first, a1, rtol=2e-10, atol=2e-10)
    assert np.allclose(second, a2, rtol=2e-7, atol=2e-7)


def test_future_branch_is_zero_from_endpoint_plus_one_without_sampling_jets() -> None:
    def unavailable_jet(j: int, x: float, y: float, z: float) -> np.ndarray:
        raise AssertionError("outside-support evaluation must not sample a jet")

    extension = BorelRightExtension(unavailable_jet, lambda j: 0, endpoint=1.0)
    assert np.array_equal(extension(0.1, -0.2, 0.3, 2.0), np.zeros(3))
    assert np.array_equal(extension(0.1, -0.2, 0.3, 4.0), np.zeros(3))


def test_glue_uses_closed_past_at_the_endpoint() -> None:
    past = lambda x, y, z, t: np.array([t, x, y + z])
    jet = lambda j, x, y, z: np.array([9.0, 8.0, 7.0]) if j == 0 else np.zeros(3)
    extension = BorelRightExtension(jet, lambda j: 1, endpoint=1.0)
    glued = extension.glue_to_past(past)

    assert np.array_equal(glued(0.2, 0.3, 0.4, 1.0), past(0.2, 0.3, 0.4, 1.0))
    assert not np.array_equal(glued(0.2, 0.3, 0.4, 1.01), past(0.2, 0.3, 0.4, 1.01))


def test_invalid_scale_and_jet_data_fail_closed() -> None:
    jet = lambda j, x, y, z: np.zeros(3)
    with pytest.raises(ValueError, match="nonnegative integer"):
        DoublingEnvelope(lambda j: 1.5)(0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        DoublingEnvelope(lambda j: -1)(0)

    bad = BorelRightExtension(
        lambda j, x, y, z: np.array([np.nan, 0.0, 0.0]),
        lambda j: 1,
    )
    with pytest.raises(ValueError, match="finite length-3 vector"):
        bad(0.0, 0.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="t>=endpoint"):
        BorelRightExtension(jet, lambda j: 1)(0.0, 0.0, 0.0, 0.9)
