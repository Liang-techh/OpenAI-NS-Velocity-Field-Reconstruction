"""Compare this checkout to a LOCAL baseline checkout; never download or alter it.

Usage: python scripts/benchmark_baseline.py /path/to/baseline --output result.json
Timings are local repeated medians, not cross-machine performance guarantees.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import timeit
import warnings
import numpy as np
from openai_ns_reconstruction.coordinates import solve_q
from openai_ns_reconstruction.profiles import LeadingProfile, toy_gaussian_profile
from openai_ns_reconstruction.status import SOURCE_PINS

EXPECTED_BLOBS={
    'coordinates.py':'28e3693804c0fa350ca4cf3ddbb5e9cdc37b166a',
    'profiles.py':'79beec1476aa668b5fba6985888eecdf041de5c5',
}


def load(path:Path,name:str):
    content=path.read_bytes()
    actual=hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()
    if actual!=EXPECTED_BLOBS[path.name]:
        raise ValueError(f'baseline file does not match the pinned snapshot: {path.name}')
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def benchmark(baseline:Path) -> dict:
    root=baseline/'src'/'openai_ns_reconstruction'
    old_c=load(root/'coordinates.py','_baseline_coordinates')
    old_p=load(root/'profiles.py','_baseline_profiles')
    root_cases=[]
    for requested_tau in (1e-6,1e-10,1e-12,1e-14):
        t=1-requested_tau
        tau=1-t
        expected=tau/(1-.6**2)
        z=expected**.495*.6
        old,new=old_c.solve_q(z,t,.005),solve_q(z,t,.005)
        root_cases.append({'tau_represented':tau,'eta':.6,'h':.005,'expected_q':expected,
            'baseline_q':old,'new_q':new,'baseline_relative_error':abs(old/expected-1),
            'new_relative_error':abs(new/expected-1)})
    p=toy_gaussian_profile()
    old_profile=old_p.LeadingProfile(p.E,p.U,p.dU_deta)
    new_profile=LeadingProfile(p.E,p.U,p.dU_deta)
    expected=.6*(-math.expm1(-2))/2
    methods={
        'baseline_trapezoid_801':lambda:old_profile.radial_average_U(2.,.6),
        'gauss_legendre_32':lambda:new_profile.radial_average_U(2.,.6),
        'optional_analytic_primitive':lambda:p.radial_average_U(2.,.6),
    }
    results={}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',DeprecationWarning)
        for name,fn in methods.items():
            fn()  # warm quadrature cache
            timings=timeit.repeat(fn,number=1000,repeat=5)
            results[name]={'value':fn(),'absolute_error':abs(fn()-expected),
                          'seconds_per_1000_calls_median':statistics.median(timings),
                          'seconds_per_1000_calls_repeats':timings}
    return {
        'kind':'local-benchmark-not-proof','sources':SOURCE_PINS,
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'platform':platform.system()},
        'near_singularity_root_cases':root_cases,
        'quadrature':{'integrand':'U(X,eta)=eta*exp(-X)','X':2.,'eta':.6,'exact_average':expected,
                      'methods':results,'number':1000,'repeat':5,'timing_statistic':'median'},
        'limitations':['One smooth integrand; no universal quadrature accuracy or timing guarantee.',
                      'Both numeric methods use the same callbacks; the optional analytic path is separate.',
                      'The benchmark does not validate the paper-specific profiles.'],
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    report=benchmark(args.baseline)
    text=json.dumps(report,indent=2,allow_nan=False)+'\n'
    if args.output:
        with args.output.open('x',encoding='utf-8') as stream:stream.write(text)
    else:print(text,end='')
