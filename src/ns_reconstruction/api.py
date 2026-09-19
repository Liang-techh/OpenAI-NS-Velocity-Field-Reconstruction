"""Chunked public API for immutable release fields and independent validation."""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import hashlib, json
import numpy as np
DATA=Path(__file__).with_name('data')

def verify_integrity()->dict:
    info=json.loads((DATA/'manifest.json').read_text())
    if info['pde_validated'] or info['paper_exact'] or info['blowup_proved'] or info['source_correspondence_verified']:
        raise ValueError('Unsupported scientific claim in manifest')
    for ident,r in info['candidates'].items():
        if ident not in ('ST054-Q2','ST054-M3') or r['pde_validated']:raise ValueError('Invalid release entry')
        for key,hkey in [('file','sha256'),('references','reference_sha256')]:
            p=(DATA/r[key]).resolve()
            if p.parent!=DATA.resolve():raise ValueError('Data path escape')
            if hashlib.sha256(p.read_bytes()).hexdigest()!=r[hkey]:raise ValueError(f'Checksum mismatch: {ident}/{key}')
    return info

class Candidate:
    """Immutable parameter storage; scalar or per-point physical times supported."""
    def __init__(self,ident:str):
        info=verify_integrity()
        if ident not in info['candidates']:raise ValueError(f'Unknown candidate: {ident}')
        from ._runtime.spacetime import Family
        self.id=ident;self.path=DATA/info['candidates'][ident]['file']
        self._family,self._raw=Family.load(self.path);self._raw.setflags(write=False)
        self.metadata=info['candidates'][ident].copy()
    @property
    def coefficients(self):
        result=self._raw.copy();result.setflags(write=False);return result
    def _inputs(self,points,t):
        x=np.asarray(points,dtype=float)
        if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():raise ValueError('Expected finite (n,3) points')
        tt=np.broadcast_to(np.asarray(t,dtype=float),(len(x),))
        if not np.isfinite(tt).all() or np.any(tt<.25) or np.any(tt>.75):raise ValueError('Time must be in [0.25,0.75]')
        return x,tt
    def fields(self,points,t):
        x,tt=self._inputs(points,t);u=np.empty_like(x);p=np.empty(len(x))
        for i in range(0,len(x),512):
            q=slice(i,i+512);u[q],p[q]=self._family.fields(self._raw,x[q],tt[q])
        return u,p
    def velocity(self,points,t):return self.fields(points,t)[0]
    def pressure(self,points,t):return self.fields(points,t)[1]
    def forcing(self,points,t):
        from ._runtime.spacetime import force
        x,tt=self._inputs(points,t);return force(x,tt,*self._raw[-2:])
    def residual(self,points,t):
        """Analytic diagnostic residual; independent acceptance uses Cartesian FD."""
        x,tt=self._inputs(points,t);out=np.empty_like(x)
        for i in range(0,len(x),256):
            q=slice(i,i+256);out[q]=self._family.analytic_residual(self._raw,x[q],tt[q])
        return out
    def validate(self,out: str|Path,seed:int=9175491)->dict:
        """Original full independent FD/quadrature protocol; no fitting occurs."""
        from ._runtime.validate import validate
        p=Path(out)
        if p.exists():raise FileExistsError('Refusing to overwrite an existing report')
        p.parent.mkdir(parents=True,exist_ok=True)
        return validate(self.path,p,seed)

@lru_cache(maxsize=2)
def load_candidate(ident:str='ST054-Q2')->Candidate:
    return Candidate(ident)
