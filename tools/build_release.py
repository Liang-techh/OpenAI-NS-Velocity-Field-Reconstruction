"""Restore a small, pinned release; never fits coefficients or edits Git refs.

Run once by the import workflow, or locally with --local-science/--local-viewer.
Every imported byte is verified BEFORE writing. No download is needed by users
of the resulting checkout. Repeated invocation refuses differing existing data.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, urllib.request
from pathlib import Path
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]

def sha(b: bytes)->str:return hashlib.sha256(b).hexdigest()
def put(path: Path, b: bytes)->None:
    if path.exists() and path.read_bytes()!=b:
        raise FileExistsError(f'Refusing to replace different release data: {path}')
    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(b)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--local-science',type=Path);ap.add_argument('--local-viewer',type=Path);args=ap.parse_args()
    manifest=json.loads((ROOT/'tools/source_manifest.json').read_text())
    receipt=[]
    for entry in manifest['entries']:
        name=entry['path'];path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT):raise ValueError('Import path escapes release root')
        if args.local_science and args.local_viewer:
            base=args.local_viewer if entry['source_path'].startswith('visualization/matlab/') else args.local_science
            rel=entry['source_path'].removeprefix('visualization/matlab/') if base==args.local_viewer else entry['source_path']
            b=(base/rel).read_bytes()
        else:
            url=f"https://raw.githubusercontent.com/{entry['source_repository']}/{entry['source_commit']}/{entry['source_path']}"
            req=urllib.request.Request(url,headers={'User-Agent':'ST054-reproducible-release'})
            with urllib.request.urlopen(req,timeout=90) as r:b=r.read()
        if sha(b)!=entry['sha256']:raise ValueError(f'Upstream checksum mismatch: {name}')
        if entry['transform']=='relative_imports':
            text=b.decode('utf-8')
            for module in ('spacetime','hybrid_basis','asymmetric_basis'):
                text=text.replace('from '+module+' import','from .'+module+' import')
            b=text.encode('utf-8')
        put(path,b)
        receipt.append(dict(path=name,upstream_sha256=entry['sha256'],release_sha256=sha(b),transform=entry['transform']))
    put(ROOT/'src/ns_reconstruction/_runtime/__init__.py',b'"""Pinned numerical implementation; see docs/PROVENANCE.md."""\n')
    sys.path.insert(0,str(ROOT/'src'))
    from ns_reconstruction._runtime.spacetime import Family
    mat_path=ROOT/'visualization/matlab/data/st054_models.mat'
    mat=loadmat(mat_path,simplify_cells=True)
    models=mat['models'];models=[models] if isinstance(models,dict) else list(models)
    candidate_meta={}
    for m in models:
        ident=str(m['id'])
        if ident not in ('ST054-Q2','ST054-M3') or bool(m['pde_validated']):raise ValueError('Unexpected model/claim')
        c=np.asarray(m['raw_coefficients'],float).ravel()
        if c.shape!=(2594,) or not np.isfinite(c).all():raise ValueError('Invalid candidate shape')
        obj=dict(schema='root_st001_compact_spacetime_v1',nr=9,nz=12,nt=8,basis_kind='hybrid9_axial_opposite3_v1',
                 coefficients=c.tolist(),metadata=dict(id=ident,source_mat_sha256=sha(mat_path.read_bytes()),
                 original_research_json_sha256=str(m['original_raw_sha256']),
                 note='Frozen MATLAB-export coefficient array. JSON metadata and small ancestor-replay rounding differ from original research JSON. No fitting or physical rescaling.'),
                 pde_validated=False,paper_exact=False,blowup_proved=False,field_identity_claim=False)
        path=ROOT/f'src/ns_reconstruction/data/{ident}.json'
        put(path,(json.dumps(obj,indent=2)+'\n').encode());f,raw=Family.load(path)
        points=np.asarray(m['ref_points']);times=np.asarray(m['ref_times']).ravel();err=np.zeros(3)
        refs=dict(points=points.tolist(),times=times.tolist(),velocity=np.asarray(m['ref_velocity']).tolist(),pressure=np.asarray(m['ref_pressure']).tolist(),residual=np.asarray(m['ref_residual']).tolist())
        for k,t in enumerate(times):
            u,p=f.fields(raw,points,t);r=f.analytic_residual(raw,points,t)
            err=np.maximum(err,[np.max(abs(u-m['ref_velocity'][k])),np.max(abs(p-m['ref_pressure'][k])),np.max(abs(r-m['ref_residual'][k]))])
        if np.max(err)>1e-8:raise ValueError(f'Imported mathematical field differs from frozen references: {ident}: {err}')
        refpath=ROOT/f'src/ns_reconstruction/data/{ident}_references.json';put(refpath,(json.dumps(refs,separators=(',',':'))+'\n').encode())
        candidate_meta[ident]=dict(file=path.name,sha256=sha(path.read_bytes()),coefficient_array_sha256=sha(c.astype('<f8').tobytes()),
                  original_research_json_sha256=str(m['original_raw_sha256']),references=refpath.name,reference_sha256=sha(refpath.read_bytes()),
                  role='primary: lower recorded volume L2' if ident=='ST054-Q2' else 'alternative: lower recorded sampled peak',
                  max_reference_errors=err.tolist(),pde_validated=False)
    catalog=dict(schema='st054_release_v1',default_candidate='ST054-Q2',candidates=candidate_meta,source_mat_sha256=sha(mat_path.read_bytes()),
                 contract=dict(nu=.01,time=[.25,.75],physical_domain='R3',evaluation_box=[[-2,2],[-2,2],[-2,2]],support='r<2 and |z|<2 for both velocity and pressure',
                 force='Original independently specified, divergence-free two-parameter curl force; a,c in [0,10]',initial_energy=1.,momentum_max_threshold=.001,momentum_volume_L2_threshold=.001,divergence_threshold=1e-5),
                 pde_validated=False,paper_exact=False,source_correspondence_verified=False,blowup_proved=False)
    put(ROOT/'src/ns_reconstruction/data/manifest.json',(json.dumps(catalog,indent=2)+'\n').encode())
    put(ROOT/'evidence/import_receipt.json',(json.dumps(dict(schema='release_import_receipt_v1',imports=receipt,candidates=candidate_meta),indent=2)+'\n').encode())
    print(json.dumps({'imported_files':len(receipt),'candidates':candidate_meta,'pde_validated':False},indent=2))
if __name__=='__main__':main()
