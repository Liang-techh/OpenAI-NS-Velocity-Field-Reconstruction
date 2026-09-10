import math
import numpy as np
import pytest

from openai_ns_reconstruction.phase import (
    TangentialBaseJet,
    backward_material_phase_defect,
    harmonic_is_single_valued,
    phase_normal,
    phase_normal_v_derivative,
    phase_value,
    primary_phase_parameters,
    projected_evolution_operator,
    reference_frame,
    shear_matrix,
    stress_cone_margins,
    tangent_frame,
)


def make_data(sigma=1):
    frame = reference_frame(1.0, np.array([-3.0, math.sqrt(91.0)]))
    params = primary_phase_parameters(
        epsilon=0.01, lambda0=frame.lambda0, u_star=1.5, L_s=80.0,
        sigma=sigma, R0=1.2, K=frame.K, g0=np.array([-3.0, math.sqrt(91.0)]),
    )
    jet = TangentialBaseJet(
        F=1.3, G=-0.4, F_R=0.2, F_Z=-0.1, F_T=0.05,
        G_R=-0.3, G_Z=0.08, G_T=-0.02,
    )
    return frame, params, jet


def test_reference_frame_and_cone_margin():
    frame, _, _ = make_data()
    assert np.linalg.norm(frame.N) == pytest.approx(1.0)
    assert frame.N @ frame.K == pytest.approx(0.0, abs=1e-15)
    assert frame.lambda0 > 0 and frame.c0 < 0
    # A target exactly opposite N is strictly inside the Eq. (7.1) cone.
    m1, m2 = stress_cone_margins(-2.0 * frame.N, frame)
    assert m1 > 0 and m2 == pytest.approx(1.0)


@pytest.mark.parametrize("sigma", [-1, 1])
def test_equation_7_2_rounding_makes_every_harmonic_single_valued(sigma):
    _, params, _ = make_data(sigma)
    assert params.k == math.ceil(params.epsilon ** -0.5)
    assert params.angular_integer != 0
    assert params.k * params.p == params.angular_integer
    for m in (-7, -1, 1, 5):
        assert harmonic_is_single_valued(params, m)


def test_equations_7_3_7_4_match_finite_difference_gradient():
    _, params, jet0 = make_data()
    R, Z, theta, v = 1.4, -0.25, 0.7, 3.0

    def jet(R_, Z_, T_):
        dR, dZ, dT = R_ - R, Z_ - Z, T_ - 0.3
        return TangentialBaseJet(
            F=jet0.F + jet0.F_R*dR + jet0.F_Z*dZ + jet0.F_T*dT,
            G=jet0.G + jet0.G_R*dR + jet0.G_Z*dZ + jet0.G_T*dT,
            F_R=jet0.F_R, F_Z=jet0.F_Z, F_T=jet0.F_T,
            G_R=jet0.G_R, G_Z=jet0.G_Z, G_T=jet0.G_T,
        )

    h = 1e-6

    def phi(R_, Z_, theta_):
        return phase_value(R_, Z_, theta_, v, params, jet(R_, Z_, 0.3))

    dR = (phi(R+h, Z, theta)-phi(R-h, Z, theta))/(2*h)
    dtheta = (phi(R, Z, theta+h)-phi(R, Z, theta-h))/(2*h)/R
    dZ = params.epsilon*(phi(R, Z+h, theta)-phi(R, Z-h, theta))/(2*h)
    n = phase_normal(R, v, params, jet0)
    assert np.allclose(n, [dR, dtheta, dZ], rtol=1e-8, atol=1e-9)


def test_backward_material_phase_defect_matches_direct_chain_rule():
    _, params, jet0 = make_data()
    R, Z, T, theta, v, b = 1.4, -0.25, 0.3, 0.7, 3.0, 0.06

    def jet(R_, Z_, T_):
        return TangentialBaseJet(
            F=jet0.F + jet0.F_R*(R_-R) + jet0.F_Z*(Z_-Z) + jet0.F_T*(T_-T),
            G=jet0.G + jet0.G_R*(R_-R) + jet0.G_Z*(Z_-Z) + jet0.G_T*(T_-T),
            F_R=jet0.F_R, F_Z=jet0.F_Z, F_T=jet0.F_T,
            G_R=jet0.G_R, G_Z=jet0.G_Z, G_T=jet0.G_T,
        )

    def phi(R_, Z_, T_, theta_, v_):
        return phase_value(R_, Z_, theta_, v_, params, jet(R_, Z_, T_))

    h = 1e-6
    dv = (phi(R, Z, T, theta, v+h)-phi(R, Z, T, theta, v-h))/(2*h)
    dT = (phi(R, Z, T+h, theta, v)-phi(R, Z, T-h, theta, v))/(2*h)
    dR = (phi(R+h, Z, T, theta, v)-phi(R-h, Z, T, theta, v))/(2*h)
    dtheta = (phi(R, Z, T, theta+h, v)-phi(R, Z, T, theta-h, v))/(2*h)
    dZ = (phi(R, Z+h, T, theta, v)-phi(R, Z-h, T, theta, v))/(2*h)
    direct = dv - params.epsilon*dT + b*dR + jet0.F*dtheta + params.epsilon*jet0.G*dZ
    exact = backward_material_phase_defect(b=b, v=v, params=params, jet=jet0)
    assert direct == pytest.approx(exact, rel=2e-8, abs=2e-9)


def test_equations_7_6_7_8_preserve_moving_incompressibility_plane():
    frame, params, jet = make_data()
    R, v = 1.4, 3.0
    n = phase_normal(R, v, params, jet)
    nv = phase_normal_v_derivative(params, jet)
    K = shear_matrix(R, jet)
    A = projected_evolution_operator(n, nv, K)
    B = tangent_frame(n, frame.c0)
    assert np.linalg.matrix_rank(B) == 2
    assert np.allclose(n @ B, 0.0, atol=1e-12)
    # Eq. (7.6) is precisely the correction making d_v(n.t)=0 for t'=A t.
    assert np.allclose(n @ A, -nv, atol=1e-12)
    for column in B.T:
        assert nv @ column + n @ (A @ column) == pytest.approx(0.0, abs=1e-12)


def test_phase_domain_checks_fail_closed():
    frame, params, jet = make_data()
    with pytest.raises(ValueError):
        phase_normal(0.0, 1.0, params, jet)
    with pytest.raises(ValueError):
        tangent_frame([1.0, 0.0, 0.0], frame.c0)
    with pytest.raises(ValueError):
        primary_phase_parameters(
            epsilon=.01, lambda0=1, u_star=1, L_s=1,
            sigma=0, R0=1, K=[1, 0], g0=[1, 0],
        )
