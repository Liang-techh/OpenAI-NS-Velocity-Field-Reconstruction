import math
import numpy as np
import pytest
from openai_ns_reconstruction.charts import DyadicChart,J_G,V_R,V_T,LAMBDA_G,T_G
from openai_ns_reconstruction.cutoffs import smooth_step,smooth_cutoff
from openai_ns_reconstruction.local_field import LocalField,LocalizedField,curl_numeric
from openai_ns_reconstruction.background import (
    BackgroundCoefficient,coefficient_vector_potential_cartesian,
)
from openai_ns_reconstruction import LeadingProfile,leading_velocity_cartesian
from openai_ns_reconstruction.verify import divergence_numeric


@pytest.mark.parametrize("ell", [8,16,40,120])
def test_dyadic_roundtrip_and_paper_covering_bounds(ell):
    chart=DyadicChart(ell)
    point=(1.2,-.3,.8)
    assert np.allclose(chart.from_physical_tau(*chart.to_physical_tau(*point)),point)
    assert 1/(T_G*chart.S_star) < chart.c_i <= 1/chart.S_star
    assert np.asarray(chart.covering_matrix()).dtype == object
    assert np.allclose(np.asarray(J_G,dtype=float)@V_R,LAMBDA_G*V_R)
    assert np.allclose(np.asarray(J_G,dtype=float)@V_T,T_G*V_T)
    assert chart.derivative_coefficients(2.)['time_slow'] == -chart.epsilon
    assert chart.normalize(1.,'pressure') == chart.Q**(1+2*chart.h)


def test_covering_is_exact_and_not_overflowed_int64():
    chart=DyadicChart(120)
    matrix=chart.covering_matrix()
    assert matrix[0,0]*matrix[1,1]-matrix[0,1]*matrix[1,0] == 14**chart.covering_index
    with pytest.raises(ValueError):
        DyadicChart(3).covering_matrix()


def test_cutoff_support_flat_edges_and_symmetry():
    assert smooth_cutoff(.5) == 1
    assert smooth_cutoff(1.) == 0
    assert smooth_cutoff(.5001) == 1
    assert smooth_cutoff(.9999) == 0
    for x in np.linspace(0,1,31):
        assert smooth_step(float(x))+smooth_step(float(1-x)) == pytest.approx(1.)
    with pytest.raises(ValueError):
        smooth_cutoff(float('nan'))


def test_curl_of_localized_potential_includes_product_rule_term():
    local=LocalField.from_axisymmetric(lambda *p:np.array([0.,0.,1.]),lambda *p:0.)
    field=LocalizedField.from_axisymmetric(local,lambda r,z,t:math.exp(-r*r-z*z))
    x,y,z,t=.3,.2,.1,.4
    c=math.exp(-x*x-y*y-z*z)
    assert np.allclose(field.velocity(x,y,z,t),[-2*y*c,2*x*c,0],atol=1e-9)
    assert abs(divergence_numeric(field.velocity,x,y,z,t)) < 1e-7


def test_cutoff_short_circuits_undefined_local_fields():
    def unavailable(*p):
        raise RuntimeError('local field must not be evaluated outside support')
    localized=LocalizedField(LocalField(unavailable,unavailable),lambda *p:0.)
    assert np.array_equal(localized.velocity(0,0,0,1.1),np.zeros(3))


def test_nonaxisymmetric_swirl_is_not_automatically_divergence_free():
    # Regression of the old overly broad documentation claim.
    local=LocalField(lambda *p:np.zeros(3),lambda x,y,z,t:x*math.hypot(x,y))
    assert divergence_numeric(local.velocity,.3,.4,0,.2) == pytest.approx(-.4,abs=1e-8)


def test_streamfunction_curl_matches_poloidal_leading_velocity():
    p=LeadingProfile(E=lambda X,e:0.,U=lambda X,e:e*(1+X),
        dU_deta=lambda X,e:1+X,average_U=lambda X,e:e*(1+X/2),
        average_dU_deta=lambda X,e:1+X/2)
    c=BackgroundCoefficient(0,p)
    A=lambda x,y,z,t:coefficient_vector_potential_cartesian(x,y,z,t,c,h=.005)
    point=(.3,.2,.1,.4)
    assert np.allclose(curl_numeric(A,*point),leading_velocity_cartesian(*point,p),
                       rtol=1e-8,atol=1e-8)
    assert np.all(np.isfinite(A(0,0,0,.3)))
    with pytest.raises(ValueError):
        BackgroundCoefficient(-1,p)
