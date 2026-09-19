"""Reproduce the original two-seed scientific rejections; never treat CI as PDE pass."""
from __future__ import annotations
import argparse,json,sys,platform
from pathlib import Path
import numpy as np, scipy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from ns_reconstruction import load_candidate

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=ROOT/'evidence/replay');a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    original=json.loads((ROOT/'evidence/upstream_ST054_results.json').read_text())
    rows=[]
    for ident in ('ST054-Q2','ST054-M3'):
        field=load_candidate(ident)
        for seed in (9175491,9175492):
            out=a.out/f'{ident}_{seed}.json'
            if out.exists():raise FileExistsError(f'Refusing existing report {out}')
            report=field.validate(out,seed)
            expected=next(r for r in original['paired_results'] if r['id']==ident and r['seed']==seed)
            fine=[r for r in report['spatial_refinement'] if r['space_step']==.005]
            value=max(r['momentum_max'] for r in fine);l2=max(r['momentum_L2'] for r in fine)
            assert len(fine)==6 and not report['all_numeric_gates_pass'] and report['pde_validated'] is False
            assert report['gates']==original['common_gates']
            assert abs(value-expected['max'])<1e-5 and abs(l2-expected['volume_L2'])<1e-5
            rows.append(dict(candidate=ident,seed=seed,candidate_sha256=report['candidate_sha256'],max=value,volume_L2=l2,failed_gates=[k for k,v in report['gates'].items() if not v],reference_tolerance=1e-5))
    summary=dict(rows=rows,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                 full_replay_completed=True,pde_validated=False,scope='Reproduction of known finite-sample rejections, not newly optimized fields or continuum validation')
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
