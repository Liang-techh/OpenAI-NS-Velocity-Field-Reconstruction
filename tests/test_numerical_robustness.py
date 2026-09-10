import math
import numpy as np
import pytest
from openai_ns_reconstruction.coordinates import (
    solve_q, solve_q_from_tau, similarity_coordinates, similarity_coordinates_from_tau,
    coordinate_derivatives, weighted_profile_derivatives,
)
from openai_ns_reconstruction.profiles import LeadingProfile, toy_gaussian_profile
from openai_ns_reconstruction.velocity import (
    blowup_probe, leading_pressure, leading_velocity_cartesian,
    leading_velocity_cylindrical_from_tau,
)
from openai_ns_reconstruction.verify import loglog_slope
from openai_ns_reconstruction.quadrature import unit_rule


@pytest.mark.parametrize('q',[1e-240,1e-80,1e-14,1e-6,1.,1e80])
@pytest.mark.parametrize('eta',[0.,.2,.6,-.99])
@pytest.mark.parametrize('h',[.005,.1,.49])
def test_scaled_root_known_solution(q,eta,h):
    tau=q*(1-eta**2)
    z=q**(.5-h)*eta
    got=solve_q_from_tau(z,tau,h)
    assert math.isclose(got,q,rel_tol=8e-12)


def test_old_near_singularity_bug():
    t=1-1e-14
    tau=1-t
    q=tau/(1-.6**2)
    assert abs(solve_q(q**.495*.6,t,.005)/q-1)<1e-12


def test_tau_api_does_not_form_t():
    p=toy_gaussian_profile()
    for tau in (1e-20,1e-80,1e-240):
        v=blowup_probe(tau,.5,p)
        assert math.isclose(v.X,.5,rel_tol=2e-14)
        assert math.isclose(v.u_theta*tau**.505,p.E(.5,0),rel_tol=2e-14)


def test_d_keeps_information_near_eta_endpoint():
    s=similarity_coordinates_from_tau(.2,1.,1e-30,.005)
    assert math.isclose(s.eta,1.,rel_tol=1e-12)
    assert s.d < 1e-20
    assert s.d>0
    assert math.isclose(s.d*s.q,s.tau,rel_tol=2e-15)


@pytest.mark.parametrize('kwargs',[
    {'z':float('nan')},{'z':float('inf')},{'tau':0},{'tau':-1},{'tau':float('inf')},
    {'h':0},{'h':.5},{'h':float('nan')},{'rtol':0},{'rtol':1},{'rtol':float('nan')},
    {'max_iter':0},{'max_iter':1.2},{'max_iter':True},
])
def test_root_invalid_input(kwargs):
    args=dict(z=.2,tau=.1,h=.005)
    args.update(kwargs)
    with pytest.raises(ValueError): solve_q_from_tau(**args)


def test_solver_reports_nonconvergence():
    with pytest.raises(RuntimeError): solve_q_from_tau(.2,.1,.005,max_iter=1)


def test_solver_reports_overflow():
    with pytest.raises(OverflowError): solve_q_from_tau(1e300,.1,.005)


@pytest.mark.parametrize('t',[1.,2.,float('nan'),float('-inf')])
def test_invalid_time(t):
    with pytest.raises(ValueError): similarity_coordinates(.2,.3,t,.005)


@pytest.mark.parametrize('r',[-1,float('nan'),float('inf')])
def test_invalid_radius(r):
    with pytest.raises(ValueError): similarity_coordinates(r,.3,.5,.005)


def test_analytic_coordinate_derivatives():
    r,z,t,h=.4,.2,.7,.005
    s=similarity_coordinates(r,z,t,h)
    deriv=coordinate_derivatives(s)
    for axis,idx in [('r',0),('z',1),('t',2)]:
        plus,minus=[r,z,t],[r,z,t]
        plus[idx]+=1e-5;minus[idx]-=1e-5
        sp,sm=similarity_coordinates(*plus,h),similarity_coordinates(*minus,h)
        for attr in ('q','eta','X') if axis!='r' else ('X',):
            got=(getattr(sp,attr)-getattr(sm,attr))/(2e-5)
            assert math.isclose(got,deriv[f'{attr}_{axis}'],rel_tol=3e-7,abs_tol=1e-8)


def test_weighted_profile_chain_rule():
    r,z,t,h,b=.4,.2,.7,.005,-.72
    s=similarity_coordinates(r,z,t,h)
    f=lambda X,e:1+X**2+3*e**2+X*e
    got=weighted_profile_derivatives(s,b,f(s.X,s.eta),2*s.X+s.eta,6*s.eta+s.X)
    def physical(r,z,t):
        ss=similarity_coordinates(r,z,t,h)
        return ss.q**b*f(ss.X,ss.eta)
    for idx,axis in enumerate(('r','z','t')):
        plus,minus=[r,z,t],[r,z,t]
        plus[idx]+=1e-5;minus[idx]-=1e-5
        numeric=(physical(*plus)-physical(*minus))/(2e-5)
        assert math.isclose(got[axis],numeric,rel_tol=3e-7)


@pytest.mark.parametrize('X',[0.,1e-240,1e-12,.1,1.,10.])
def test_gaussian_radial_average(X):
    p=toy_gaussian_profile()
    numeric=LeadingProfile(p.E,p.U,p.dU_deta)
    assert math.isclose(numeric.radial_average_U(X,.6),p.radial_average_U(X,.6),rel_tol=2e-14)
    assert math.isclose(numeric.radial_average_dU_deta(X,.6),p.radial_average_dU_deta(X,.6),rel_tol=2e-14)


def test_quadrature_polynomial_and_readonly_cache():
    p=LeadingProfile(lambda X,e:0.,lambda X,e:X**13+e,lambda X,e:1.)
    assert math.isclose(p.radial_average_U(1.3,.2,n=8),1.3**13/14+.2,rel_tol=5e-14)
    x,w=unit_rule(8)
    with pytest.raises(ValueError): x[0]=1
    assert unit_rule(8)[0] is x


@pytest.mark.parametrize('n',[0,1,2.5,True,4096])
def test_invalid_quadrature_order(n):
    with pytest.raises(ValueError): unit_rule(n)


def test_regular_pressure_on_axis_and_radial_balance():
    p=toy_gaussian_profile()
    assert p.pressure_radial_derivative(0,.3)==p.F(0,.3)**2
    X,e=.4,.3
    diff=(p.Pi(X+1e-6,e)-p.Pi(X-1e-6,e))/2e-6
    assert math.isclose(diff,p.pressure_radial_derivative(X,e),rel_tol=1e-9)
    r,z,t=.4,.2,.7
    pressure_r=(leading_pressure(r+1e-6,0,z,t,p)-leading_pressure(r-1e-6,0,z,t,p))/2e-6
    u=leading_velocity_cartesian(r,0,z,t,p)
    assert math.isclose(pressure_r,u[1]**2/r,rel_tol=1e-8)


def test_axis_limit_requires_profile_information():
    p=LeadingProfile(lambda X,e:0.,lambda X,e:0.,lambda X,e:0.)
    with pytest.raises(ValueError): p.pressure_radial_derivative(0,0)
    with pytest.raises(ValueError): leading_pressure(.1,0,0,.5,p)


def test_nonzero_small_radius_is_not_forced_to_zero():
    v=leading_velocity_cylindrical_from_tau(1e-16,0,.1,toy_gaussian_profile())
    assert v.u_r!=0
    assert v.u_theta!=0


def test_invalid_axis_profile_not_silently_hidden():
    p=LeadingProfile(lambda X,e:1.,lambda X,e:0.,lambda X,e:0.)
    with pytest.raises(ValueError): leading_velocity_cartesian(0,0,0,.5,p)


def test_profile_cannot_claim_exact_without_provenance():
    with pytest.raises(ValueError):
        LeadingProfile(lambda X,e:0.,lambda X,e:0.,lambda X,e:0.,paper_exact=True)


@pytest.mark.parametrize('xs,ys',[
    ([],[]),([1],[2]),([1,1],[2,3]),([1,2],[1]),([1,2],[1,float('nan')]),
    ([0,1],[1,2]),([1,2],[0,2]),([[1,2]],[[3,4]]),([1,float('inf')],[2,3]),
])
def test_bad_loglog_samples(xs,ys):
    with pytest.raises(ValueError): loglog_slope(np.array(xs),np.array(ys))
