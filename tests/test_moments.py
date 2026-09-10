import math
import numpy as np
import pytest
from openai_ns_reconstruction.moments import solve_quadratic_moments, power_moment_matrix


def test_scalar_quadratic_known_small_root():
    # c+c^2=d has a stable closed form 2d/(1+sqrt(1+4d)).
    d=.05
    out=solve_quadratic_moments(np.eye(1),np.ones((1,1,1)),np.array([d]))
    expected=2*d/(1+math.sqrt(1+4*d))
    assert math.isclose(out.coefficients[0],expected,rel_tol=1e-12)
    assert out.residual_norm<1e-13
    assert out.smallness_estimate==.4
    assert out.contraction_estimate<=.5
    assert out.iterations>1
    with pytest.raises(ValueError): out.coefficients[0]=0


def test_linear_system_has_no_smallness_requirement():
    B=np.array([[2.,1.],[0.,3.]])
    d=np.array([1e10,-1e10])
    out=solve_quadratic_moments(B,np.zeros((2,2,2)),d)
    assert out.iterations==1
    assert np.allclose(B@out.coefficients,d)
    assert out.smallness_estimate==0


def test_nonsymmetric_bilinear_system_known_coefficients():
    B=np.array([[2.,.2],[.1,3.]])
    Q=np.array([[[.1,.3],[.2,-.1]],[[0.,.2],[.1,.15]]])
    expected=np.array([.01,-.015])
    d=B@expected+np.einsum('ijk,j,k->i',Q,expected,expected)
    out=solve_quadratic_moments(B,Q,d)
    assert np.allclose(out.coefficients,expected,rtol=1e-12,atol=1e-15)
    assert np.linalg.norm(out.coefficients)<=out.ball_radius_estimate


def test_zero_rhs_and_tiny_rhs():
    for d in (0.,1e-100,1e-240):
        out=solve_quadratic_moments(np.eye(1),np.ones((1,1,1)),np.array([d]))
        assert math.isclose(out.coefficients[0],d,rel_tol=1e-12,abs_tol=0)


def test_insufficient_smallness_is_not_silent_success():
    with pytest.raises(ValueError,match='smallness'):
        solve_quadratic_moments(np.eye(1),np.ones((1,1,1)),np.array([1.]))


@pytest.mark.parametrize('B',[np.zeros((2,2)),np.array([[1.,1.],[1.,1.]])])
def test_singular_matrix_rejected(B):
    with pytest.raises(ValueError,match='singular'):
        solve_quadratic_moments(B,np.zeros((2,2,2)),np.ones(2))


def test_max_iterations_respected():
    with pytest.raises(RuntimeError):
        solve_quadratic_moments(np.eye(1),np.ones((1,1,1)),np.array([.05]),max_iter=1)


@pytest.mark.parametrize('kwargs',[{'rtol':0},{'atol':-1},{'max_iter':0},{'max_iter':1.5},
                                 {'max_iter':True},{'rtol':float('nan')}])
def test_bad_controls(kwargs):
    with pytest.raises(ValueError):
        solve_quadratic_moments(np.eye(1),np.ones((1,1,1)),np.array([.05]),**kwargs)


def test_bad_shapes_and_nonfinite():
    with pytest.raises(ValueError): solve_quadratic_moments(np.eye(2),np.zeros((1,1,1)),np.ones(2))
    with pytest.raises(ValueError): solve_quadratic_moments(np.eye(1),np.zeros((1,1,1)),np.array([np.nan]))


def test_power_moment_matrix_normalization_and_symmetry():
    intervals=np.array([[1.,2.],[3.,4.]])
    B=power_moment_matrix(np.array([0.,1.]),intervals)
    assert np.allclose(B,np.array([[1.,1.],[1.5,3.5]]),atol=1e-14)
    assert np.linalg.det(B)>0
    d=np.array([.1,.2])
    out=solve_quadratic_moments(B,np.zeros((2,2,2)),d)
    assert np.allclose(B@out.coefficients,d,atol=1e-14)


@pytest.mark.parametrize('exps,intervals',[
    ([0.,0.],[[1,2],[3,4]]),([0,1],[[1,3],[2,4]]),([0,1],[[-1,0],[1,2]]),
    ([0,1],[[3,4],[1,2]]),([],[]),([float('nan')],[[1,2]]),
])
def test_bad_moment_geometry(exps,intervals):
    with pytest.raises(ValueError): power_moment_matrix(exps,intervals)
