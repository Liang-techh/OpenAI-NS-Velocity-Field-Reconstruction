"""Reproducible diagnostic report; no test here certifies a complete paper field."""
from __future__ import annotations
import math
import platform
import numpy as np
from .status import construction_status
from .coordinates import solve_q_from_tau
from .profiles import toy_gaussian_profile
from .velocity import blowup_probe
from .background import (BackgroundCoefficient, background_potential_cartesian,
    background_velocity_cartesian, background_swirl_cartesian)
from .local_field import curl_numeric
from .moments import solve_quadratic_moments
from .verify import divergence_numeric, navier_stokes_residual_numeric, loglog_slope


def run_diagnostics() -> dict:
    profile=toy_gaussian_profile()
    root_errors=[]
    for q in (1e-240,1e-80,1e-14,1e-6,1.,1e80):
        for eta in (0.,.2,.6,-.99):
            z=q**.495*eta
            root_errors.append(abs(solve_q_from_tau(z,q*(1-eta*eta),.005)/q-1))
    taus=np.logspace(-2,-80,40)
    velocities=np.array([blowup_probe(float(tau),.5,profile).u_theta for tau in taus])
    slope=loglog_slope(taus,velocities)
    coefficients=[BackgroundCoefficient(0,profile),BackgroundCoefficient(2,profile,2.5)]
    point=(.25,-.17,.12,.7)
    potential=lambda x,y,z,t:background_potential_cartesian(x,y,z,t,coefficients)
    velocity=lambda x,y,z,t:background_velocity_cartesian(x,y,z,t,coefficients)
    numeric=curl_numeric(potential,*point,eps=1e-5)+background_swirl_cartesian(*point,coefficients)
    curl_error=float(np.linalg.norm(velocity(*point)-numeric))
    divergence=abs(divergence_numeric(velocity,*point,eps=1e-5))
    nu=.37
    def u(x,y,z,t):
        return math.exp(-2*nu*t)*np.array([math.sin(x)*math.cos(y),-math.cos(x)*math.sin(y),0.])
    def p(x,y,z,t):
        return math.exp(-4*nu*t)*(math.cos(2*x)+math.cos(2*y))/4
    steps=(.05,.025,.0125)
    residuals=[float(np.linalg.norm(navier_stokes_residual_numeric(u,p,.7,.4,.2,.3,
                viscosity=nu,eps_space=e,eps_time=e))) for e in steps]
    ratios=[residuals[i]/residuals[i+1] for i in range(2)]
    moment=solve_quadratic_moments(np.eye(1),np.ones((1,1,1)),np.array([.05]))
    checks={
        "similarity_max_relative_error": {"value":max(root_errors),"tolerance":1e-12},
        "toy_exponent_error": {"value":abs(slope+.505),"tolerance":1e-12},
        "background_curl_error": {"value":curl_error,"tolerance":1e-7},
        "background_divergence": {"value":divergence,"tolerance":1e-7},
        "taylor_green_finest_residual": {"value":residuals[-1],"tolerance":1e-4},
        "taylor_green_refinement_ratio_error": {"value":max(abs(v-4) for v in ratios),"tolerance":.2},
        "pointwise_moment_residual": {"value":moment.residual_norm,"tolerance":1e-12},
    }
    for check in checks.values():
        check["passed"]=bool(math.isfinite(check["value"]) and check["value"]<=check["tolerance"])
    return {
        "schema_version":1,"kind":"numerical-diagnostics-not-proof",
        "diagnostics_passed":all(c["passed"] for c in checks.values()),
        "environment":{"python":platform.python_version(),"numpy":np.__version__,"platform":platform.system()},
        "parameters":{"h":.005,"profile":profile.name,"quadrature_nodes":32,
                      "root_grid_cases":len(root_errors),"toy_probe_samples":len(taus),
                      "background_orders":[0,2],"background_cutoff_scale":2.5,
                      "background_point":point,"background_eps_space":1e-5,
                      "taylor_green_point":[.7,.4,.2,.3],"viscosity":nu},
        "checks":checks,"measured_toy_exponent":slope,
        "refinement":{"steps":steps,"residual_norms":residuals,"ratios":ratios},
        "construction":construction_status(),
    }
