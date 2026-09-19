from __future__ import annotations
import argparse,json
from . import load_candidate,verify_integrity

def main():
    p=argparse.ArgumentParser(description='Frozen ST054 fields; full NS acceptance remains unmet')
    p.add_argument('command',choices=['status','verify','evaluate','validate'])
    p.add_argument('--candidate',default='ST054-Q2',choices=['ST054-Q2','ST054-M3'])
    p.add_argument('--point',nargs=3,type=float,default=[.1,0,.1]);p.add_argument('--time',type=float,default=.5)
    p.add_argument('--out');p.add_argument('--seed',type=int,default=9175491)
    a=p.parse_args()
    if a.command in ('status','verify'):
        print(json.dumps(verify_integrity(),indent=2));return 0
    f=load_candidate(a.candidate)
    if a.command=='evaluate':
        u,pres=f.fields([a.point],a.time);r=f.residual([a.point],a.time)
        print(json.dumps(dict(candidate=a.candidate,time=a.time,point=a.point,velocity=u[0].tolist(),pressure=float(pres[0]),analytic_residual=r[0].tolist(),pde_validated=False),indent=2));return 0
    if not a.out:p.error('--out is required for validation')
    r=f.validate(a.out,a.seed);failed=[k for k,v in r['gates'].items() if not v]
    print('FAIL: '+', '.join(failed) if failed else 'Sampled numerical gates pass; not a continuum certificate')
    return int(bool(failed))
if __name__=='__main__':raise SystemExit(main())
