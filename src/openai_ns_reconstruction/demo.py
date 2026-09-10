"""Reproducible diagnostic artifacts. Gaussian samples are explicitly TOY data."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
import scipy
from .coordinates import similarity_coordinates_from_tau
from .profiles import LeadingProfile
from .velocity import blowup_probe, leading_velocity_cartesian
from .heat_exterior import HeatExterior
from .provenance import status_report
from .verify import loglog_slope, navier_stokes_residual_numeric


def toy_gaussian_profile() -> LeadingProfile:
    def radial_average(X):
        return 1.0 if X == 0 else -math.expm1(-X) / X
    return LeadingProfile(
        E=lambda X,e:math.sqrt(2*X)*math.exp(-X)*(1+.1*e*e),
        U=lambda X,e:e*math.exp(-X), dU_deta=lambda X,e:math.exp(-X),
        Pi=lambda X,e:-.5*math.exp(-2*X)*(1+.1*e*e)**2,
        name="toy-gaussian-not-openai-profile", paper_exact=False,
        average_U=lambda X,e:e*radial_average(X),
        average_dU_deta=lambda X,e:radial_average(X),
        F=lambda X,e:math.exp(-X)*(1+.1*e*e))


def write_json(path: Path, value: dict) -> None:
    # Reports are small; serialize fully before opening so NaN cannot leave a partial report.
    text = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def generate_demo(output: Path) -> dict:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    h, X_in = .005, .5
    profile = toy_gaussian_profile()
    taus = np.logspace(-2,-30,15)
    values = [blowup_probe(float(tau),X_in,profile,h=h).u_theta for tau in taus]
    with (output/'toy_blowup_probe.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['tau','u_theta','profile_status'])
        writer.writerows((float(tau),float(v),'toy') for tau,v in zip(taus,values))
    with (output/'toy_velocity_samples.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['x','y','z','t','u_x','u_y','u_z','profile_status'])
        for t in (.5,.9,.99):
            for z in (-.03,0.,.03):
                for x in np.linspace(-.2,.2,5):
                    for y in np.linspace(-.2,.2,5):
                        v=leading_velocity_cartesian(float(x),float(y),z,t,profile,h=h)
                        writer.writerow([float(x),float(y),z,t,*v.tolist(),'toy'])
    heat=HeatExterior(h=h,c_inf=1.)
    with (output/'heat_exterior.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['radius','tau','K','scope'])
        for tau in (1.,.1,.01,0.):
            for r in (.5,1.,2.,4.):
                writer.writerow([r,tau,heat.swirl_from_tau(r,tau),'paper-formula-exterior-only'])
    coordinate_errors=[]
    for q in (1.,1e-12,1e-30,1e-100,1e-200):
        for eta in (-.9,.2,.99):
            s=similarity_coordinates_from_tau(math.sqrt(q),q**(.5-h)*eta,q*(1-eta**2),h)
            coordinate_errors.append(abs(s.q/q-1))
    nu=.7
    def tg_u(x,y,z,t):
        return math.exp(-2*nu*t)*np.array([math.sin(x)*math.cos(y),-math.cos(x)*math.sin(y),0.])
    def tg_p(x,y,z,t):
        return math.exp(-4*nu*t)*(math.cos(2*x)+math.cos(2*y))/4
    steps=(.04,.02,.01)
    tg_errors=[float(np.linalg.norm(navier_stokes_residual_numeric(
        tg_u,tg_p,.7,.4,.1,.3,viscosity=nu,eps_space=e,eps_time=e))) for e in steps]
    heat_defect=float(np.linalg.norm(navier_stokes_residual_numeric(
        heat.velocity,heat.pressure,1.4,.3,0.,.6,eps_space=2e-4,eps_time=2e-4)))
    metrics={"max_coordinate_q_relative_error":max(coordinate_errors),
             "toy_loglog_slope":loglog_slope(taus,np.array(values)),
             "toy_expected_slope":-.5-h,
             "heat_ode_defect_max":max(abs(heat.ode_defect(Z)) for Z in (0.,.1,1.,10.,100.)),
             "heat_exterior_ns_residual_norm":heat_defect,
             "independent_taylor_green_steps":list(steps),
             "independent_taylor_green_residual_norms":tg_errors}
    checks={"coordinate_accuracy":metrics['max_coordinate_q_relative_error']<2e-12,
            "toy_scaling":abs(metrics['toy_loglog_slope']+.5+h)<1e-12,
            "heat_ode":metrics['heat_ode_defect_max']<1e-10,
            "heat_exterior_balance":heat_defect<2e-6,
            "independent_refinement":tg_errors[1]<tg_errors[0]/3.8 and tg_errors[2]<tg_errors[1]/3.8}
    report=status_report()
    report.update({"artifact_kind":"numerical-diagnostics-not-full-reconstruction",
                   "parameters":{"h":h,"X_in":X_in,"heat_c_inf":1.},
                   "profile_status":"toy", "metrics":metrics,"checks":checks,
                   "all_diagnostic_checks_passed":all(checks.values()),
                   "environment":{"python":platform.python_version(),"numpy":np.__version__,
                                  "scipy":scipy.__version__},
                   "artifacts":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(output.glob('*.csv'))},
                   "package_source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                             for p in sorted(Path(__file__).parent.glob('*.py'))},
                   "limitations":["Gaussian samples are not the paper's constructed leading profile.",
                       "Heat exterior is valid away from the axis only and is not the global field.",
                       "No smooth extension of the full residual through t=1 is certified.",
                       "No Lean compilation or full-paper theorem cross-check was performed."]})
    write_json(output/'report.json',report)
    return report
