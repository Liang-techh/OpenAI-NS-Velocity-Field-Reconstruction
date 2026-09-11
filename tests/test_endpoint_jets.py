import numpy as np
import pytest

from openai_ns_reconstruction.endpoint_jets import (
    borel_extension_from_full_spacetime_jets,
    normal_jet_family,
    normal_time_jet_from_full,
    validate_full_spacetime_jet,
)


def _manufactured_full_jets(degree: int, x: float, y: float, z: float) -> np.ndarray:
    """Endpoint t=1 jets of F=(t+x,t^2+y,t^3+z)."""
    shape = (3,) + (4,) * degree
    tensor = np.zeros(shape)
    if degree == 0:
        tensor[:] = [1.0 + x, 1.0 + y, 1.0 + z]
    elif degree == 1:
        tensor[:, 0] = [1.0, 2.0, 3.0]
        tensor[0, 1] = 1.0
        tensor[1, 2] = 1.0
        tensor[2, 3] = 1.0
    elif degree == 2:
        tensor[:, 0, 0] = [0.0, 2.0, 6.0]
    elif degree == 3:
        tensor[:, 0, 0, 0] = [0.0, 0.0, 6.0]
    return tensor


def test_full_jet_normal_contraction_selects_repeated_time_direction() -> None:
    point = (0.2, -0.3, 0.4)
    expected = (
        np.array([1.2, 0.7, 1.4]),
        np.array([1.0, 2.0, 3.0]),
        np.array([0.0, 2.0, 6.0]),
        np.array([0.0, 0.0, 6.0]),
    )
    for degree, want in enumerate(expected):
        got = normal_time_jet_from_full(_manufactured_full_jets, degree, *point)
        assert np.array_equal(got, want)


def test_normal_family_is_directly_compatible_with_borel_extension() -> None:
    point = (0.2, -0.3, 0.4)
    normal = normal_jet_family(_manufactured_full_jets)
    assert np.array_equal(normal(2, *point), np.array([0.0, 2.0, 6.0]))

    extension = borel_extension_from_full_spacetime_jets(
        _manufactured_full_jets, lambda degree: 1, endpoint=1.0
    )
    t = 1.04
    expected = np.array([t + point[0], t * t + point[1], t**3 + point[2]])

    # The doubling scales are 1,2,4,8,16,... .  At s=0.04 the degree 0--3
    # terms all lie on the exact cutoff plateau; degree 4 is identically zero
    # for this polynomial and degree 5 is outside support before its jet is
    # sampled.  Thus this checks the adapter and Borel composition against the
    # analytic polynomial, rather than reconstructing the production sum.
    assert np.allclose(extension(*point, t), expected, rtol=0.0, atol=2e-15)


def test_spatial_derivative_entries_are_not_mistaken_for_normal_time_jets() -> None:
    tensor = np.zeros((3, 4))
    tensor[:, 0] = [2.0, 3.0, 4.0]
    tensor[:, 1:] = [[100.0, 200.0, 300.0], [400.0, 500.0, 600.0], [700.0, 800.0, 900.0]]
    full = lambda degree, x, y, z: tensor
    assert np.array_equal(
        normal_time_jet_from_full(full, 1, 0.0, 0.0, 0.0),
        np.array([2.0, 3.0, 4.0]),
    )


def test_full_jet_shape_and_finiteness_fail_closed() -> None:
    with pytest.raises(ValueError, match="shape"):
        validate_full_spacetime_jet(np.zeros((3, 4)), 2)
    bad = np.zeros((3, 4, 4))
    bad[1, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        validate_full_spacetime_jet(bad, 2)
    with pytest.raises(ValueError, match="nonnegative integer"):
        validate_full_spacetime_jet(np.zeros(3), True)
    with pytest.raises(ValueError, match="nonnegative integer"):
        validate_full_spacetime_jet(np.zeros(3), -1)
    with pytest.raises(ValueError, match="callable"):
        normal_jet_family(None)
