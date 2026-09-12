"""Status, fail-closed audit, reproducible verification, and labelled toy exports."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
from .status import construction_status
from .diagnostics import run_diagnostics
from .profiles import toy_gaussian_profile
from .velocity import leading_velocity_cartesian_from_tau


def construction_audit() -> dict:
    """Return a compact audit view derived only from the runtime truth surface."""
    status = construction_status()
    return {
        "schema_version": 1,
        "kind": "construction-audit",
        "paper_exact_velocity_available": status["paper_exact_velocity_available"],
        "construction_status": status["status"],
        "completion_gate_passed": status["paper_exact_velocity_available"],
        "sources": status["sources"],
        "stages": [
            {
                "id": stage["id"],
                "name": stage["name"],
                "status": stage["status"],
                "remaining": stage["remaining"],
            }
            for stage in status["stages"]
        ],
        "limitations": status["limitations"],
    }


def export_demo(output: Path, *, grid_size: int=41, tau: float=1e-3, h: float=.005) -> dict:
    if isinstance(grid_size,bool) or not isinstance(grid_size,int) or not 5<=grid_size<=401:
        raise ValueError("grid_size must be an integer between 5 and 401")
    if not math.isfinite(tau) or tau<=0 or not math.isfinite(h) or not 0<h<.01:
        raise ValueError("demo requires tau>0 and 0<h<1/100, all finite")
    output=Path(output)
    if output.exists():
        raise FileExistsError("output directory exists; choose a new directory")
    profile=toy_gaussian_profile()
    x=np.linspace(-3*math.sqrt(tau),3*math.sqrt(tau),grid_size)
    z=np.linspace(-.8*tau**(.5-h),.8*tau**(.5-h),grid_size)
    velocity=np.empty((grid_size,grid_size,3))
    for j,zj in enumerate(z):
        for i,xi in enumerate(x):
            velocity[j,i]=leading_velocity_cartesian_from_tau(float(xi),0.,float(zj),tau,profile,h=h)
    output.mkdir(parents=True,exist_ok=False)
    data=output/'toy_velocity_slice.npz'
    np.savez_compressed(data,x=x,z=z,velocity=velocity,tau=tau,h=h,
                        truth_status='toy-not-openai',profile=profile.name)
    manifest={
        "schema_version":1,"kind":"toy-not-openai",
        "paper_exact_velocity_available":False,
        "description":"Cartesian components on a y=0 slice; velocity[j,i,k] matches z[j],x[i],component k.",
        "parameters":{"tau":tau,"h":h,"grid_size":grid_size,"profile":profile.name},
        "artifact":{"name":data.name,"sha256":hashlib.sha256(data.read_bytes()).hexdigest(),
                    "velocity_shape":list(velocity.shape)},
        "construction":construction_status(),
    }
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return manifest


def main(argv: list[str] | None=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    status=sub.add_parser('status',help='Print construction status (not a completion claim).')
    status.add_argument('--require-paper-exact',action='store_true')
    audit=sub.add_parser('audit',help='Print a fail-closed construction/provenance audit view.')
    audit.add_argument('--require-paper-exact',action='store_true')
    verify=sub.add_parser('verify',help='Run independent numerical diagnostics.')
    verify.add_argument('--output',type=Path,help='Write a new JSON report; existing files are not overwritten.')
    verify.add_argument('--require-paper-exact',action='store_true')
    demo=sub.add_parser('demo',help='Export a labelled toy slice, NOT the paper field.')
    demo.add_argument('--output',required=True,type=Path)
    demo.add_argument('--grid-size',type=int,default=41)
    demo.add_argument('--tau',type=float,default=1e-3)
    demo.add_argument('--h',type=float,default=.005)
    demo.add_argument('--require-paper-exact',action='store_true')
    args=parser.parse_args(argv)
    try:
        if (
            args.command == 'demo'
            and args.require_paper_exact
            and not construction_status()['paper_exact_velocity_available']
        ):
            print('Incomplete paper construction: no paper-exact velocity instance is available.',file=sys.stderr)
            return 2
        if args.command=='status':
            result=construction_status()
        elif args.command=='audit':
            result=construction_audit()
        elif args.command=='verify':
            result=run_diagnostics()
            if args.output:
                with args.output.open('x',encoding='utf-8') as stream:
                    stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
        else:
            result=export_demo(args.output,grid_size=args.grid_size,tau=args.tau,h=args.h)
        print(json.dumps(result,indent=2,allow_nan=False))
        if getattr(args,'require_paper_exact',False) and not construction_status()['paper_exact_velocity_available']:
            print('Incomplete paper construction: no paper-exact velocity instance is available.',file=sys.stderr)
            return 2
        if args.command=='verify' and not result['diagnostics_passed']:
            return 1
        return 0
    except (ValueError,RuntimeError,OSError,OverflowError) as exc:
        parser.exit(2,f'error: {exc}\n')


if __name__=='__main__':
    raise SystemExit(main())
