import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from openai_ns_reconstruction.cli import main, export_demo
from openai_ns_reconstruction.diagnostics import run_diagnostics
from openai_ns_reconstruction.status import construction_status


def test_status_does_not_promote_toy_or_tests_to_complete(capsys):
    assert main(['status'])==0
    status=json.loads(capsys.readouterr().out)
    assert status['paper_exact_velocity_available'] is False
    assert len(status['stages'])==9
    status['stages'][0]['status']='fake'
    assert construction_status()['stages'][0]['status']!='fake'


def test_runtime_status_tracks_current_formal_structure_without_promotion():
    status=construction_status()
    stages={stage['id']:stage for stage in status['stages']}
    assert status['paper_exact_velocity_available'] is False
    assert status['sources']['pr5_integration_merge_commit']=='c0a08e68b0f65a59a4dca70c50934957a0608475'
    assert status['sources']['lean_source_reviewed'] is True
    assert status['sources']['lean_compiled'] is False
    assert 'NaturalAxisRange' in stages[1]['implemented']
    assert 'Eq. (5.5)' in stages[2]['implemented']
    assert 'PhaseEstimates' in stages[3]['implemented']
    assert 'analytic cutoff gradient' in stages[7]['implemented']
    assert stages[1]['status']==stages[2]['status']==stages[3]['status']==stages[7]['status']=='formal-structure'
    assert stages[8]['status']=='diagnostic-only'


@pytest.mark.parametrize('command',['status','verify'])
def test_strict_gate_fails_closed(command,capsys):
    assert main([command,'--require-paper-exact'])==2
    output=capsys.readouterr()
    assert 'Incomplete paper construction' in output.err
    result=json.loads(output.out)
    if command=='verify':
        assert result['diagnostics_passed']
        assert not result['construction']['paper_exact_velocity_available']


def test_verify_writes_reproducible_data_not_a_certificate(tmp_path,capsys):
    path=tmp_path/'verification.json'
    assert main(['verify','--output',str(path)])==0
    result=json.loads(path.read_text())
    assert result['kind']=='numerical-diagnostics-not-proof'
    assert result['diagnostics_passed']
    assert result['environment']['numpy']==np.__version__
    assert result==json.loads(capsys.readouterr().out)
    with pytest.raises(SystemExit) as exc:
        main(['verify','--output',str(path)])
    assert exc.value.code==2


def test_demo_export_label_shape_and_checksum(tmp_path):
    path=tmp_path/'demo'
    manifest=export_demo(path,grid_size=9,tau=1e-30)
    data=path/'toy_velocity_slice.npz'
    assert manifest['kind']=='toy-not-openai'
    assert hashlib.sha256(data.read_bytes()).hexdigest()==manifest['artifact']['sha256']
    with np.load(data,allow_pickle=False) as values:
        assert values['velocity'].shape==(9,9,3)
        assert np.all(np.isfinite(values['velocity']))
        assert str(values['truth_status'])=='toy-not-openai'
        assert values['x'][4]==0
    with pytest.raises(FileExistsError): export_demo(path)


@pytest.mark.parametrize('kwargs',[{'tau':0},{'tau':float('nan')},{'h':.02},{'h':float('inf')},
                                 {'grid_size':3},{'grid_size':True},{'grid_size':1000}])
def test_invalid_demo_requests_do_not_write(tmp_path,kwargs):
    target=tmp_path/'invalid'
    with pytest.raises(ValueError): export_demo(target,**kwargs)
    assert not target.exists()


def test_module_entrypoint_subprocess():
    import os
    env=dict(os.environ,PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
    result=subprocess.run([sys.executable,'-m','openai_ns_reconstruction','status'],
                          env=env,capture_output=True,text=True,check=True)
    assert not json.loads(result.stdout)['paper_exact_velocity_available']


def test_diagnostics_report_is_deterministic():
    assert run_diagnostics()==run_diagnostics()
