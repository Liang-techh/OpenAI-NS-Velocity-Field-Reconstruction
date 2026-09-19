"""Publication regression tests; passing is not PDE acceptance."""
import json,sys
from pathlib import Path
import numpy as np
import pytest
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ns_candidate import load,verify

@pytest.fixture(scope='module')
def models():return loadmat(ROOT/'visualization/matlab/data/st054_models.mat',simplify_cells=True)['models']

def test_integrity_and_english():
    verify()
    import re
    for path in ROOT.rglob('*'):
        if path.suffix in ('.md','.py','.m','.json','.yml','.txt') and not any(p in path.parts for p in ('.git','outputs','__pycache__')):
            assert not re.search('[\u3400-\u9fff]',path.read_text(encoding='utf-8')),str(path)

@pytest.mark.parametrize('ident',['ST054-Q2','ST054-M3'])
def test_native_reference_parity(models,ident):
    model=next(m for m in models if m['id']==ident);f,raw=load(ident)
    np.testing.assert_array_equal(raw,model['raw_coefficients'])
    for i,t in enumerate(model['ref_times']):
        u,p=f.fields(raw,model['ref_points'],t);r=f.analytic_residual(raw,model['ref_points'],t)
        np.testing.assert_allclose(u,model['ref_velocity'][i],atol=1e-8,rtol=1e-8)
        np.testing.assert_allclose(p,model['ref_pressure'][i],atol=1e-8,rtol=1e-8)
        np.testing.assert_allclose(r,model['ref_residual'][i],atol=1e-8,rtol=1e-8)

@pytest.mark.parametrize('ident',['ST054-Q2','ST054-M3'])
def test_independent_residual_and_support(ident):
    from spacetime import force,quad
    from validate import cartesian_residual
    f,raw=load(ident);xyz=np.array([[.1,.2,.3],[.2,0,-.2],[0,0,.1],[.3,0,1.7]])
    r,d=cartesian_residual(lambda x,t:f.fields(raw,x,t),lambda x,t:force(x,t,*raw[-2:]),xyz,.517,.00125,.000625)
    np.testing.assert_allclose(r,f.analytic_residual(raw,xyz,.517),atol=1e-6,rtol=1e-4)
    assert abs(d).max()<1e-7
    for a in f.fields(raw,np.array([[2,0,0],[0,0,-2],[0,0,2],[3,3,3]]),.5):assert np.count_nonzero(a)==0
    s,z,w=quad(96);u,_=f.fields(raw,np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    with pytest.raises(ValueError):f.fields(raw,xyz,.751)

def test_report_gates():
    reports=list((ROOT/'evidence/replay').glob('*.json'));assert len(reports)==4
    for p in reports:
        r=json.loads(p.read_text());assert not r['pde_validated'] and not r['all_numeric_gates_pass']
        assert not r['gates']['momentum_max'] and not r['gates']['momentum_L2']
