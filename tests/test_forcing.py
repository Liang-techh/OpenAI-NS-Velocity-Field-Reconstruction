import numpy as np

from openai_ns_reconstruction.verify import (
    divergence_numeric,
    forced_ns_closure_error_numeric,
    reconstruct_forcing_numeric,
)


def solid_rotation(x: float, y: float, z: float, t: float) -> np.ndarray:
    """Divergence-free solid-body rotation u=(-y,x,0)."""
    return np.array([-y, x, 0.0])


def zero_pressure(x: float, y: float, z: float, t: float) -> float:
    return 0.0


def exact_force(x: float, y: float, z: float, t: float) -> np.ndarray:
    # For u=(-y,x,0): u_t=0, Delta u=0, (u.grad)u=(-x,-y,0).
    return np.array([-x, -y, 0.0])


def test_solid_rotation_is_divergence_free() -> None:
    div = divergence_numeric(solid_rotation, 0.7, -0.4, 0.2, 0.3, eps=1e-6)
    assert abs(div) < 1e-9


def test_reconstructed_force_matches_manufactured_solution() -> None:
    point = (0.7, -0.4, 0.2, 0.3)
    got = reconstruct_forcing_numeric(
        solid_rotation,
        zero_pressure,
        *point,
        viscosity=0.37,
        eps_space=1e-5,
        eps_time=1e-6,
    )
    expected = exact_force(*point)
    assert np.allclose(got, expected, rtol=1e-7, atol=1e-7)


def test_forced_ns_closure_error_vanishes_for_exact_force() -> None:
    point = (-0.3, 0.8, -0.1, 0.6)
    defect = forced_ns_closure_error_numeric(
        solid_rotation,
        zero_pressure,
        exact_force,
        *point,
        viscosity=1.2,
        eps_space=1e-5,
        eps_time=1e-6,
    )
    assert np.linalg.norm(defect) < 1e-7
