import math
import numpy as np
import pytest
from openai_ns_reconstruction.heat_exterior import HeatExterior, heat_factor
from openai_ns_reconstruction.coordinates import similarity_coordinates_from_tau
from openai_ns_reconstruction.verify import navier_stokes_residual_numeric


@pytest.mark.parametrize("m", range(6))
def test_heat_endpoint_derivative_formula(m):
    h=.005
    exact=(-1)**m*math.prod(h+k for k in range(m))*math.prod(1+h+k for k in range(m))
    assert heat_factor(0,h,derivative=m) == exact


@pytest.mark.parametrize("Z", [0.,1e-6,.1,1.,10.,100.])
def test_heat_factor_positivity_and_ode(Z):
    H=heat_factor(Z)
    assert 0 < H <= 1
    H1=heat_factor(Z,derivative=1)
    assert H1 < 0
    assert 0 <= -Z*H1/H < .005
    assert abs(HeatExterior().ode_defect(Z)) < 1e-10


def test_similarity_heat_profile_matches_physical_swirl_and_endpoints():
    heat=HeatExterior()
    for r,z,tau in [(1.,.4,.2),(2.,-.8,.01),(.1,.05,1e-8)]:
        s=similarity_coordinates_from_tau(r,z,tau,heat.h)
        profile=s.q**(-s.A)*heat.profile_E(s.X,s.eta)
        assert profile == pytest.approx(heat.swirl_from_tau(r,tau),rel=1e-11)
    assert heat.profile_E(5.,1.) == pytest.approx(5**(-heat.A))
    assert heat.swirl_from_tau(2.,0) == pytest.approx(2**(-heat.A))


def test_heat_swirl_pde_by_independent_finite_differences():
    heat=HeatExterior()
    r,t=1.3,.6
    errors=[]
    for e in [.02,.01,.005]:
        K=heat.swirl(r,t)
        kt=(heat.swirl(r,t+e)-heat.swirl(r,t-e))/(2*e)
        kr=(heat.swirl(r+e,t)-heat.swirl(r-e,t))/(2*e)
        krr=(heat.swirl(r+e,t)-2*K+heat.swirl(r-e,t))/e**2
        errors.append(abs(kt-krr-kr/r+K/r**2))
    assert errors[-1] < errors[0]/12
    assert errors[-1] < 1e-4


def test_exterior_pressure_and_full_ns_balance():
    heat=HeatExterior()
    r,t,e=1.4,.6,1e-4
    dp=(heat.pressure_from_tau(r+e,1-t)-heat.pressure_from_tau(r-e,1-t))/(2*e)
    assert dp == pytest.approx(heat.swirl(r,t)**2/r, rel=2e-7)
    defect=navier_stokes_residual_numeric(heat.velocity,heat.pressure,
        1.4,.3,0.,.6,eps_space=2e-4,eps_time=2e-4)
    assert np.linalg.norm(defect) < 2e-6


def test_heat_exterior_cannot_masquerade_as_regular_axis_profile():
    with pytest.raises(ValueError):
        HeatExterior().swirl(0,.5)
    with pytest.raises(ValueError):
        heat_factor(-1)
    with pytest.raises(ValueError):
        HeatExterior().profile_E(1,1.1)
