"""Finite-difference diagnostics, not proofs of smoothness or NS blow-up.

R(u,p)=u_t+(u.grad)u-nu*Delta u+grad p is the force required by a field.
R(u,p)-f is an independent closure test ONLY when f is independently supplied.
In particular, defining f by this same routine is a tautology, not verification.
"""
from __future__ import annotations
from typing import Callable
import math
import numpy as np
from .coordinates import _finite
from .quadrature import unit_rule

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]
TimeDomain = tuple[float | None, float | None] | None


def finite_vector(value: object) -> np.ndarray:
    v = np.asarray(value,dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)):
        raise ValueError("field must return a finite length-3 vector")
    return v


def _step(center: float, eps: float) -> tuple[float,float]:
    center, eps = _finite(center,"coordinate"), _finite(eps,"eps")
    if eps <= 0:
        raise ValueError("finite-difference step must be positive")
    minus, plus = center-eps, center+eps
    if not (math.isfinite(minus) and math.isfinite(plus) and minus<center<plus):
        raise ValueError("finite-difference stencil is not representable")
    return minus, plus


def _sample(u: VectorField, point: tuple[float,...] | list[float]) -> np.ndarray:
    if not all(math.isfinite(float(v)) for v in point):
        raise ValueError("sample point must be finite")
    return finite_vector(u(*point))


def jacobian_numeric(u: VectorField, x: float, y: float, z: float, t: float, *,
                     eps: float = 1e-5) -> np.ndarray:
    J = np.zeros((3,3))
    for axis in range(3):
        minus,plus = _step((x,y,z)[axis],eps)
        p,m = [x,y,z,t],[x,y,z,t]
        p[axis],m[axis] = plus,minus
        J[:,axis] = (_sample(u,p)-_sample(u,m))/(plus-minus)
    return J


def divergence_numeric(u: VectorField, x: float, y: float, z: float, t: float, *,
                       eps: float = 1e-5) -> float:
    return float(np.trace(jacobian_numeric(u,x,y,z,t,eps=eps)))


def laplacian_vector_numeric(u: VectorField, x: float, y: float, z: float, t: float, *,
                            eps: float = 1e-4) -> np.ndarray:
    center = _sample(u,(x,y,z,t))
    out = np.zeros(3)
    for axis in range(3):
        c = (x,y,z)[axis]
        minus,plus = _step(c,eps)
        p,m = [x,y,z,t],[x,y,z,t]
        p[axis],m[axis] = plus,minus
        hp,hm = plus-c,c-minus
        # Also valid for the slightly unequal offsets created by floating point.
        out += 2*((_sample(u,p)-center)/hp-(center-_sample(u,m))/hm)/(hp+hm)
    return finite_vector(out)


def time_derivative_numeric(u: VectorField, x: float, y: float, z: float, t: float, *,
                            eps: float = 1e-5, time_domain: TimeDomain = None) -> np.ndarray:
    t,eps = _finite(t,"t"),_finite(eps,"eps")
    if eps <= 0:
        raise ValueError("eps must be positive")
    if time_domain is not None:
        if len(time_domain) != 2:
            raise ValueError("time_domain is (lower, upper), with optional None endpoints")
        lower,upper = time_domain
        if lower is not None:
            lower = _finite(lower,"lower time bound")
            if t <= lower:
                raise ValueError("center must be inside the open time domain")
            eps = min(eps,(t-lower)/4)
        if upper is not None:
            upper = _finite(upper,"upper time bound")
            if t >= upper:
                raise ValueError("center must be inside the open time domain")
            eps = min(eps,(upper-t)/4)
    minus,plus = _step(t,eps)
    if time_domain is not None:
        lower,upper = time_domain
        if (lower is not None and minus<=lower) or (upper is not None and plus>=upper):
            raise ValueError("time stencil crosses the field domain; use analytic derivatives")
    return (_sample(u,(x,y,z,plus))-_sample(u,(x,y,z,minus)))/(plus-minus)


def gradient_scalar_numeric(p: ScalarField, x: float, y: float, z: float, t: float, *,
                            eps: float = 1e-5) -> np.ndarray:
    _finite(t,"t")
    out = np.zeros(3)
    for axis in range(3):
        minus,plus = _step((x,y,z)[axis],eps)
        xp,xm = [x,y,z,t],[x,y,z,t]
        xp[axis],xm[axis] = plus,minus
        out[axis] = (_finite(p(*xp),"pressure")-_finite(p(*xm),"pressure"))/(plus-minus)
    return out


def navier_stokes_terms_numeric(u: VectorField, p: ScalarField,
    x: float, y: float, z: float, t: float, *, viscosity: float = 1.0,
    eps_space: float = 1e-4, eps_time: float = 1e-5,
    time_domain: TimeDomain = None) -> dict[str,np.ndarray]:
    """Expose each momentum term so cancellation cannot hide a sign mistake."""
    viscosity = _finite(viscosity,"viscosity")
    if viscosity < 0:
        raise ValueError("viscosity must be nonnegative")
    uv = _sample(u,(x,y,z,t))
    return {"time":time_derivative_numeric(u,x,y,z,t,eps=eps_time,time_domain=time_domain),
            "advection":jacobian_numeric(u,x,y,z,t,eps=eps_space)@uv,
            "viscous":-viscosity*laplacian_vector_numeric(u,x,y,z,t,eps=eps_space),
            "pressure":gradient_scalar_numeric(p,x,y,z,t,eps=eps_space)}


def navier_stokes_residual_numeric(u: VectorField, p: ScalarField,
    x: float, y: float, z: float, t: float, *, viscosity: float = 1.0,
    eps_space: float = 1e-4, eps_time: float = 1e-5,
    time_domain: TimeDomain = None) -> np.ndarray:
    terms = navier_stokes_terms_numeric(u,p,x,y,z,t,viscosity=viscosity,
        eps_space=eps_space,eps_time=eps_time,time_domain=time_domain)
    return finite_vector(sum(terms.values(),np.zeros(3)))


def reconstruct_forcing_numeric(u: VectorField, p: ScalarField,
    x: float, y: float, z: float, t: float, *, viscosity: float = 1.0,
    eps_space: float = 1e-4, eps_time: float = 1e-5,
    time_domain: TimeDomain = None) -> np.ndarray:
    return navier_stokes_residual_numeric(u,p,x,y,z,t,viscosity=viscosity,
        eps_space=eps_space,eps_time=eps_time,time_domain=time_domain)


def forced_ns_closure_error_numeric(u: VectorField, p: ScalarField, forcing: VectorField,
    x: float, y: float, z: float, t: float, *, viscosity: float = 1.0,
    eps_space: float = 1e-4, eps_time: float = 1e-5,
    time_domain: TimeDomain = None) -> np.ndarray:
    return navier_stokes_residual_numeric(u,p,x,y,z,t,viscosity=viscosity,
        eps_space=eps_space,eps_time=eps_time,time_domain=time_domain)-_sample(forcing,(x,y,z,t))


def loglog_slope(xs: np.ndarray, ys: np.ndarray) -> float:
    xs,ys = np.asarray(xs,dtype=float),np.asarray(ys,dtype=float)
    if xs.ndim!=1 or xs.shape!=ys.shape or len(xs)<2:
        raise ValueError("xs and ys must have equal 1D shapes and at least two samples")
    if not np.all(np.isfinite(xs)) or not np.all(np.isfinite(ys)) or np.any(xs<=0) or np.any(ys==0):
        raise ValueError("finite xs>0 and nonzero finite ys are required")
    lx,ly = np.log(xs),np.log(np.abs(ys))
    lx,ly = lx-lx.mean(),ly-ly.mean()
    denominator = float(lx@lx)
    if denominator <= 0:
        raise ValueError("xs must contain at least two distinct values")
    return float(lx@ly/denominator)


def kinetic_energy_axisymmetric(u: Callable[[float,float,float],np.ndarray], t: float, *,
    radius: float, z_min: float, z_max: float, n: int = 32) -> float:
    """Integral of |u|^2/2 over a FINITE cylinder, including the 2*pi*r measure.

    The three components must be cylindrical and independent of angle.
    This is neither total whole-space energy nor a uniform-in-time bound.
    """
    radius,z_min,z_max = _finite(radius,"radius"),_finite(z_min,"z_min"),_finite(z_max,"z_max")
    t = _finite(t,"t")
    if radius<=0 or z_max<=z_min:
        raise ValueError("positive radius and z_min<z_max are required")
    nodes,weights = unit_rule(n)
    total = 0.0
    for R,wR in zip(nodes,weights):
        r = radius*float(R)
        for Z,wZ in zip(nodes,weights):
            z = z_min+(z_max-z_min)*float(Z)
            v = finite_vector(u(r,z,t))
            total += float(wR*wZ)*r*float(v@v)
    return math.pi*radius*(z_max-z_min)*total
