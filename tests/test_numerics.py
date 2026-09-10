import math
import numpy as np
import pytest
from openai_ns_reconstruction import (
    LeadingProfile, solve_q, solve_q_from_tau, similarity_coordinates_from_tau,
    leading_velocity_cylindrical, leading_velocity_cartesian, leading_pressure_cartesian,
)
from openai_ns_reconstruction.velocity import blowup_probe
from openai_ns_reconstruction.verify import (
    divergence_numeric, loglog_slope, time_derivative_numeric,
    navier_stokes_residual_numeric, gradient_scalar_numeric,
)


@pytest.mark.parametrize("q", [1.0, 1e-12, 1e-30, 1e-100, 1e-200])
@pytest.mark.parametrize("eta", [-0.9, 0.0, 0.2, 0.99])
def test_small_scale_coordinates_recover_manufactured_root(q, eta):
    h, X = .005, .7
    tau = q * (1 - eta**2)
    s = similarity_coordinates_from_tau(math.sqrt(q) * math.sqrt(2 * X),
                                       q**(.5 - h) * eta, tau, h)
    assert s.q == pytest.approx(q, rel=2e-12, abs=0)
    assert s.X == pytest.approx(X, rel=2e-12)
    assert s.eta == pytest.approx(eta, rel=2e-12, abs=1e-15)
    assert s.d == pytest.approx(1 - eta**2, rel=2e-12, abs=0)


@pytest.mark.parametrize("args", [(float('nan'),1,.005), (0,0,.005), (0,1,0),
                                  (0,1,.5), (0,float('inf'),.005)])
def test_bad_coordinate_inputs(args):
    with pytest.raises(ValueError):
        solve_q_from_tau(*args)


def test_time_api_cannot_hide_lost_tau():
    with pytest.raises(ValueError):
        solve_q(.1, 1.0, .005)
    assert solve_q_from_tau(0, 1e-40, .005) == 1e-40
    with pytest.raises(ValueError):
        solve_q_from_tau(.1, .1, .005, rtol=0)
    with pytest.raises(RuntimeError):
        solve_q_from_tau(.4, .1, .005, max_iter=1)


def polynomial_profile(**kwargs):
    return LeadingProfile(E=lambda X,e: math.sqrt(2*X), U=lambda X,e:e*X**4,
                          dU_deta=lambda X,e:X**4, F=lambda X,e:1., **kwargs)


def test_gauss_average_and_axis_limits():
    p = polynomial_profile()
    for X in [0, 1e-20, .5, 3.]:
        assert p.radial_average_U(X,.7,n=8) == pytest.approx(.7*X**4/5, rel=4e-15, abs=1e-15)
        assert p.radial_average_dU_deta(X,.7,n=8) == pytest.approx(X**4/5, rel=4e-15, abs=1e-15)
    assert p.pressure_radial_derivative(0,.1) == 1
    with pytest.raises(ValueError):
        p.radial_average_U(1,.2,n=1)


def test_exact_average_callbacks_skip_quadrature():
    p = polynomial_profile(average_U=lambda X,e:e*X**4/5,
                           average_dU_deta=lambda X,e:X**4/5)
    assert p.radial_average_U(3,.4) == pytest.approx(.4*81/5)


def test_nonzero_axis_swirl_is_rejected_not_silently_erased():
    bad = LeadingProfile(lambda X,e:1., lambda X,e:0., lambda X,e:0.)
    with pytest.raises(ValueError):
        leading_velocity_cylindrical(0,0,.5,bad)
    with pytest.raises(ValueError):
        leading_velocity_cylindrical(-.1,0,.5,polynomial_profile())
    with pytest.raises(ValueError):
        leading_pressure_cartesian(.1,.2,0,.5,polynomial_profile())


def test_very_small_nonzero_radius_is_not_artificially_zeroed():
    v = leading_velocity_cylindrical(1e-16,0,.5,polynomial_profile())
    assert v.u_theta > 0


def test_probe_preserves_tau_far_below_time_resolution():
    p = polynomial_profile()
    taus = np.logspace(-20,-100,9)
    vals = [blowup_probe(float(t), .7, p).u_theta for t in taus]
    assert loglog_slope(taus, vals) == pytest.approx(-.505, abs=1e-12)


def test_leading_divergence_with_polynomial_averages():
    p = polynomial_profile()
    u = lambda x,y,z,t: leading_velocity_cartesian(x,y,z,t,p)
    assert abs(divergence_numeric(u,.2,.1,.07,.3,eps=1e-5)) < 1e-7


def test_stencils_stay_inside_time_domain():
    sampled = []
    def u(x,y,z,t):
        assert 0 <= t < 1
        sampled.append(t)
        return np.array([t*t, t, 1.])
    assert np.allclose(time_derivative_numeric(u,0,0,0,0,eps=.01,time_bounds=(0,1)),
                       [0,1,0], atol=1e-12)
    t = 1-1e-9
    value = time_derivative_numeric(u,0,0,0,t,eps=.01,time_bounds=(0,1))
    assert np.allclose(value,[2*t,1,0],rtol=1e-6)
    assert min(sampled) >= 0 and max(sampled) < 1


def taylor_green(viscosity=.7):
    def u(x,y,z,t):
        a = math.exp(-2*viscosity*t)
        return a*np.array([math.sin(x)*math.cos(y), -math.cos(x)*math.sin(y),0])
    def p(x,y,z,t):
        return math.exp(-4*viscosity*t)*(math.cos(2*x)+math.cos(2*y))/4
    return u,p


def test_independent_viscous_nonlinear_pressure_balance_and_refinement():
    # This solution has nonzero u_t, advection, viscosity and pressure gradient;
    # its force is analytically zero, NOT computed from the residual routine.
    u,p = taylor_green()
    errors = [np.linalg.norm(navier_stokes_residual_numeric(
        u,p,.7,.4,.1,.3,viscosity=.7,eps_space=e,eps_time=e)) for e in (.04,.02,.01)]
    assert errors[1] < errors[0]/3.8
    assert errors[2] < errors[1]/3.8
    assert errors[-1] < 1e-4


@pytest.mark.parametrize("xs,ys", [([1,1],[1,2]), ([1],[2]), ([1,2],[1,np.nan]),
                                    ([1,2],[0,1]), ([1,2],[1,2,3])])
def test_slope_rejects_degenerate_data(xs,ys):
    with pytest.raises(ValueError):
        loglog_slope(xs,ys)


def test_bad_stencils_and_vector_shapes():
    with pytest.raises(ValueError):
        divergence_numeric(lambda *p:[1,2],0,0,0,.2)
    with pytest.raises(ValueError):
        divergence_numeric(lambda *p:[1,2,3],0,0,0,.2,eps=0)
    with pytest.raises(ValueError):
        gradient_scalar_numeric(lambda *p:float('nan'),0,0,0,.2)
