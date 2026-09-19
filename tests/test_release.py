from pathlib import Path
import json
import numpy as np
import pytest
from ns_reconstruction import load_candidate,verify_integrity
ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('ident',['ST054-Q2','ST054-M3'])
def test_references(ident):
    f=load_candidate(ident);d=verify_integrity()['candidates'][ident]
    refs=json.loads((f.path.parent/d['references']).read_text())
    for k,t in enumerate(refs['times']):
        u,p=f.fields(refs['points'],t);r=f.residual(refs['points'],t)
        np.testing.assert_allclose(u,refs['velocity'][k],atol=1e-8,rtol=0)
        np.testing.assert_allclose(p,refs['pressure'][k],atol=1e-8,rtol=0)
        np.testing.assert_allclose(r,refs['residual'][k],atol=1e-8,rtol=0)

@pytest.mark.parametrize('ident',['ST054-Q2','ST054-M3'])
def test_physical_contract(ident):
    f=load_candidate(ident);r=f.coefficients
    assert not r.flags.writeable and len(r)==2594
    with pytest.raises(ValueError):r[0]=0
    x=[[2,0,0],[0,0,2],[0,0,-2],[3,0,0]]
    for v in f.fields(x,.5):assert np.max(abs(v))==0
    from ns_reconstruction._runtime.spacetime import quad
    s,z,w=quad(96);u=f.velocity(np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    pts=np.array([[.23,.15,.37],[0,0,-.19]])
    a=.7;Q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    np.testing.assert_allclose(f.velocity(pts@Q.T,.51),f.velocity(pts,.51)@Q.T,atol=2e-10)

@pytest.mark.parametrize('ident',['ST054-Q2','ST054-M3'])
def test_independent_operator(ident):
    from ns_reconstruction._runtime.validate import cartesian_residual
    f=load_candidate(ident);x=np.array([[.17,.13,.18],[.3,-.1,-.33],[0,0,.1]])
    r,div=cartesian_residual(f.fields,f.forcing,x,.537,.00125,.000625)
    np.testing.assert_allclose(r,f.residual(x,.537),atol=2e-6,rtol=0)
    assert np.max(abs(div))<1e-7

def test_invalid_and_empty_inputs():
    f=load_candidate()
    for t in [-1,.751,float('nan')]:
        with pytest.raises(ValueError):f.fields([[0,0,0]],t)
    with pytest.raises(ValueError):f.fields([[np.nan,0,0]],.5)
    with pytest.raises(ValueError):load_candidate('ST055')
    assert f.fields(np.empty((0,3)),.5)[0].shape==(0,3)
    t=np.array([.25,.5,.75]);x=np.tile([.1,0,.1],(3,1));u,_=f.fields(x,t)
    assert not np.allclose(u[0],u[-1])

def test_matlab_data_identity():
    from scipy.io import loadmat
    models=loadmat(ROOT/'visualization/matlab/data/st054_models.mat',simplify_cells=True)['models']
    for m in models:
        np.testing.assert_array_equal(load_candidate(str(m['id'])).coefficients,m['raw_coefficients'])

def test_all_text_is_english_and_local_links_exist():
    import re
    for p in ROOT.rglob('*'):
        if p.is_file() and p.suffix in ('.md','.py','.m','.json','.yml','.txt') and '.git' not in p.parts and 'outputs' not in p.parts:
            text=p.read_text(encoding='utf-8')
            assert not re.search('[\u3400-\u9fff]',text),str(p)
    for p in ROOT.rglob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in target or target.startswith('#'):continue
            assert (p.parent/target.split('#')[0]).exists(),(str(p),target)
