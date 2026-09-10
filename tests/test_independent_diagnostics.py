import math
import numpy as np
import pytest
from openai_ns_reconstruction.verify import (
    divergence_numeric, forced_ns_closure_error_numeric, navier_stokes_terms_numeric,
    navier_stokes_residual_numeric, time_derivative_numeric, kinetic_energy_axisymmetric,
)
from openai_ns_reconstruction.local_field import curl_numeric


def manufactured_u(x,y,z,t):
    return np.array([t*y*y,t*z*z,t*x*x])


def manufactured_p(x,y,z,t):
    return t*(x*x+y*y+z*z)


def manufactured_force(x,y,z,t,nu):
    return (np.array([y*y,z*z,x*x]) + 2*t*t*np.array([y*z*z,z*x*x,x*y*y])
            -2*nu*t*np.ones(3)+2*t*np.array([x,y,z]))


def test_full_independent_forced_solution_all_four_terms():
    point=(.7,-.4,.2,.3);nu=.37
    force=lambda x,y,z,t:manufactured_force(x,y,z,t,nu)
    defect=forced_ns_closure_error_numeric(manufactured_u,manufactured_p,force,*point,
        viscosity=nu,eps_space=1e-3,eps_time=1e-4)
    assert np.linalg.norm(defect)<1e-9
    terms=navier_stokes_terms_numeric(manufactured_u,manufactured_p,*point,viscosity=nu,eps_space=1e-3)
    assert all(np.linalg.norm(v)>0 for v in terms.values())
    assert abs(divergence_numeric(manufactured_u,*point))<1e-12
    wrong=lambda x,y,z,t:force(x,y,z,t)+np.array([1.,0.,0.])
    assert np.linalg.norm(forced_ns_closure_error_numeric(manufactured_u,manufactured_p,wrong,*point,
        viscosity=nu,eps_space=1e-3))>.9


def test_taylor_green_unforced_closure_and_second_order_convergence():
    nu=.37
    def u(x,y,z,t):
        a=math.exp(-2*nu*t)
        return a*np.array([math.sin(x)*math.cos(y),-math.cos(x)*math.sin(y),0.])
    def p(x,y,z,t):
        return math.exp(-4*nu*t)*(math.cos(2*x)+math.cos(2*y))/4
    errors=[]
    for eps in (.05,.025,.0125):
        errors.append(np.linalg.norm(navier_stokes_residual_numeric(u,p,.7,.4,.2,.3,
            viscosity=nu,eps_space=eps,eps_time=eps)))
    assert errors[2]<1e-4
    assert 3.8<errors[0]/errors[1]<4.2
    assert 3.8<errors[1]/errors[2]<4.2


def test_time_domain_never_crosses_singularity():
    calls=[]
    def u(x,y,z,t):
        assert t<1
        calls.append(t)
        return np.array([t,2*t,3*t])
    out=time_derivative_numeric(u,0,0,0,1-1e-10,eps=1e-3,time_domain=(None,1.))
    assert np.allclose(out,[1,2,3],atol=1e-5)
    assert len(calls)==2


@pytest.mark.parametrize('eps',[0.,-1.,float('nan'),float('inf')])
def test_invalid_steps(eps):
    with pytest.raises(ValueError): curl_numeric(manufactured_u,0,0,0,.5,eps=eps)
    with pytest.raises(ValueError): time_derivative_numeric(manufactured_u,0,0,0,.5,eps=eps)


def test_domain_and_viscosity_validation():
    with pytest.raises(ValueError): time_derivative_numeric(manufactured_u,0,0,0,1.,time_domain=(None,1.))
    with pytest.raises(ValueError): navier_stokes_residual_numeric(manufactured_u,manufactured_p,0,0,0,.5,viscosity=-1)
    with pytest.raises(ValueError): time_derivative_numeric(manufactured_u,0,0,0,1e50,eps=1e-8)


@pytest.mark.parametrize('bad',[np.zeros((3,1)),np.zeros(2),np.array([1.,float('nan'),0.])])
def test_bad_field_shape_and_nonfinite(bad):
    with pytest.raises(ValueError): divergence_numeric(lambda *x:bad,.2,.3,.4,.5)


def test_cylindrical_energy_measure():
    R,L=2.,3.
    constant=kinetic_energy_axisymmetric(lambda r,z,t:np.array([1.,2.,3.]),.3,radius=R,z_min=0,z_max=L)
    rotation=kinetic_energy_axisymmetric(lambda r,z,t:np.array([0.,r,0.]),.3,radius=R,z_min=0,z_max=L)
    assert math.isclose(constant,7*math.pi*R**2*L,rel_tol=1e-13)
    assert math.isclose(rotation,math.pi*R**4*L/4,rel_tol=1e-13)
