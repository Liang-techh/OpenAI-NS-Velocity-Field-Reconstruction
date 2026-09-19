"""Frozen-candidate CLI. Scientific validation failure is exit 1, not success."""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'runtime'))
from spacetime import Family,force

def verify():
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for path,expected in manifest['files'].items():
        p=(ROOT/path).resolve()
        if not p.is_relative_to(ROOT) or not p.is_file():raise ValueError('Missing or escaping path: '+path)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise ValueError('Checksum mismatch: '+path)
    return manifest

def load(candidate='ST054-Q2'):
    reg=json.loads((ROOT/'candidates/registry.json').read_text())
    if candidate not in reg['candidates']:raise ValueError('Unknown candidate: '+candidate)
    row=reg['candidates'][candidate];p=ROOT/row['file']
    if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Candidate hash mismatch')
    return Family.load(p)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['status','verify','evaluate','validate'])
    p.add_argument('--candidate',default='ST054-Q2',choices=['ST054-Q2','ST054-M3']);p.add_argument('--point',nargs=3,type=float,default=[.1,0,.1]);p.add_argument('--time',type=float,default=.5)
    p.add_argument('--seed',type=int,default=9175491);p.add_argument('--out',type=Path,default=ROOT/'outputs/validation.json');a=p.parse_args()
    if a.command=='status':print((ROOT/'candidates/registry.json').read_text());return 0
    if a.command=='verify':m=verify();print(json.dumps({'integrity':'passed','files':len(m['files']),'pde_validated':False},indent=2));return 0
    f,raw=load(a.candidate)
    if a.command=='evaluate':
        xyz=np.array([a.point]);u,pressure=f.fields(raw,xyz,a.time);R=f.analytic_residual(raw,xyz,a.time)
        print(json.dumps({'candidate':a.candidate,'time':a.time,'point':a.point,'velocity':u[0].tolist(),'pressure':float(pressure[0]),'force':force(xyz,a.time,*raw[-2:])[0].tolist(),'residual':R[0].tolist(),'pde_validated':False},indent=2));return 0
    if a.out.exists():p.error('Output exists; refusing to overwrite evidence')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    from validate import validate
    rep=validate(ROOT/f'candidates/{a.candidate}/candidate.json',a.out,a.seed)
    failed=[k for k,v in rep['gates'].items() if not v]
    print('FAIL: '+', '.join(failed) if failed else 'Sampled numerical gates passed; no continuum certificate')
    return int(bool(failed))
if __name__=='__main__':raise SystemExit(main())
