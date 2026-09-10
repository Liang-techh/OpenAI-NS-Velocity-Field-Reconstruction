import math
import numpy as np
import pytest
from openai_ns_reconstruction.background import (
    BackgroundCoefficient, background_velocity_cartesian, background_potential_cartesian,
    background_swirl_cartesian, background_pressure,
)
from openai_ns_reconstruction.coordinates import similarity_coordinates
from openai_ns_reconstruction.cutoffs import standard_cutoff, standard_cutoff_derivative, CompactSpatialCutoff
from openai_ns_reconstruction.local_field import LocalField, LocalizedField, curl_numeric
from openai_ns_reconstruction.profiles import LeadingProfile, toy_gaussian_profile
from openai_ns_reconstruction.velocity import leading_velocity_cartesian, leading_pressure
from openai_ns_reconstruction.verify import divergence_numeric


@pytest.mark.parametrize('s',[-1.,0.,.5,1.,2.])
def test_cutoff_outside_transition(s):
    assert standard_cutoff(s)==(1. if s<=.5 else 0.)
    assert standard_cutoff_derivative(s)==0.


@pytest.mark.parametrize('s',[.51,.6,.7,.75,.8,.9,.99])
def test_cutoff_analytic_derivative(s):
    num=(standard_cutoff(s+1e-6)-standard_cutoff(s-1e-6))/2e-6
    assert math.isclose(num,standard_cutoff_derivative(s),rel_tol=1e-6,abs_tol=1e-9)


def test_cutoff_endpoints_flat_and_range():
    for s in np.linspace(.5,1,301):
        assert 0<=standard_cutoff(s)<=1
        assert standard_cutoff_derivative(s)<=0
    for delta in (1e-3,1e-4,1e-6):
        assert abs(standard_cutoff_derivative(.5+delta))<1e-100
        assert abs(standard_cutoff_derivative(1-delta))<1e-100


def test_zeroth_background_is_leading_field():
    p=toy_gaussian_profile();coeff=[BackgroundCoefficient(0,p)]
    for point in ((.2,.3,.1,.4),(0,0,.2,.4)):
        assert np.allclose(background_velocity_cartesian(*point,coeff),leading_velocity_cartesian(*point,p),atol=1e-13)
        assert math.isclose(background_pressure(*point,coeff),leading_pressure(*point,p),rel_tol=1e-14)


@pytest.mark.parametrize('n',[1,2,5])
def test_positive_order_analytic_curl_in_cutoff_transition(n):
    p=toy_gaussian_profile();point=(.2,.3,.1,.4)
    q=similarity_coordinates(math.hypot(*point[:2]),point[2],point[3],.005).q
    coeff=[BackgroundCoefficient(n,p,.75/q)]
    A=lambda *x:background_potential_cartesian(*x,coeff)
    exact=background_velocity_cartesian(*point,coeff)
    numeric=curl_numeric(A,*point,eps=1e-5)+background_swirl_cartesian(*point,coeff)
    assert np.allclose(exact,numeric,rtol=5e-8,atol=1e-8)
    u=lambda *x:background_velocity_cartesian(*x,coeff)
    assert abs(divergence_numeric(u,*point,eps=2e-5))<1e-7


def test_coefficient_axis_and_validation():
    p=toy_gaussian_profile()
    for n in (-1,True,1.5):
        with pytest.raises(ValueError): BackgroundCoefficient(n,p)
    with pytest.raises(ValueError): BackgroundCoefficient(1,p,0)
    coeff=[BackgroundCoefficient(2,p)]
    u=background_velocity_cartesian(0,0,.2,.4,coeff)
    assert u[0]==u[1]==0
    assert np.all(np.isfinite(u))
    with pytest.raises(ValueError):
        background_velocity_cartesian(.2,.3,.1,.4,coeff,cutoff=lambda s:1.)


def test_missing_pressure_rejected():
    p=LeadingProfile(lambda X,e:0.,lambda X,e:0.,lambda X,e:0.)
    with pytest.raises(ValueError): background_pressure(.2,.3,.1,.4,[BackgroundCoefficient(0,p)])


def test_localization_product_rule_and_divergence():
    A=lambda x,y,z,t:np.array([0.,0.,-.5*(x*x+y*y)])
    curl=lambda x,y,z,t:np.array([-y,x,0.])
    cutoff=CompactSpatialCutoff(1.)
    local=LocalField(A,lambda *x:0.,analytic_curl=curl)
    field=LocalizedField(local,cutoff,cutoff.gradient)
    point=(.6,.5,.1,.3)
    got=field.velocity(*point)
    expected=curl_numeric(field.localized_potential,*point,eps=1e-6)
    assert np.allclose(got,expected,rtol=1e-8,atol=1e-8)
    assert not np.allclose(got,cutoff(*point)*curl(*point))  # omitted product term would fail
    assert abs(divergence_numeric(field.velocity,*point,eps=1e-5))<2e-7
    assert np.array_equal(field.velocity(2,0,0,.3),np.zeros(3))


def test_generic_nonaxisymmetric_cutoff_is_not_incompressible():
    # Regression against falsely claiming every Cartesian callback preserves divergence.
    local=LocalField(lambda *x:np.zeros(3),lambda x,y,z,t:math.hypot(x,y),analytic_curl=lambda *x:np.zeros(3))
    field=LocalizedField(local,lambda x,y,z,t:x,lambda *x:np.array([1.,0.,0.]))
    assert math.isclose(divergence_numeric(field.velocity,.4,.7,.2,.3),-.7,rel_tol=1e-9)
