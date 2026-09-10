import math

from openai_ns_reconstruction.coordinates import coordinate_identity_error, similarity_coordinates


def test_similarity_identities():
    h = 0.005
    for z, t in [(0.0, 0.3), (0.1, 0.7), (-0.05, 0.95), (0.003, 0.999)]:
        ez, etau = coordinate_identity_error(0.2, z, t, h)
        assert ez < 1e-10
        assert etau < 1e-10


def test_center_plane_has_q_equal_tau():
    h = 0.005
    t = 0.93
    tau = 1.0 - t
    s = similarity_coordinates(0.4, 0.0, t, h)
    assert math.isclose(s.q, tau, rel_tol=0.0, abs_tol=1e-14)
    assert s.eta == 0.0


def test_singular_path_keeps_X_fixed():
    h = 0.005
    X_in = 0.7
    for tau in (1e-2, 1e-4, 1e-6):
        r = math.sqrt(2.0 * X_in * tau)
        s = similarity_coordinates(r, 0.0, 1.0 - tau, h)
        assert math.isclose(s.X, X_in, rel_tol=1e-10, abs_tol=1e-10)
