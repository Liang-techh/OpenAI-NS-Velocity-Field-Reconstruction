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
    assert 'PressureDatum.Admissible' in stages[1]['implemented']
    assert 'S=32' in stages[1]['implemented']
    assert 'finalAngular' in stages[1]['implemented']
    assert 'clockWeight' in stages[1]['implemented']
    assert 'actual all-real-line SchedulePressure.axisPressure' in stages[1]['implemented']
    assert 'low-|Z| uniform H^2 margin' in stages[1]['implemented']
    assert 'sigma=sqrt(m)/20' in stages[1]['implemented']
    assert 'Lambda=max(1+B+L,1+(14000/9)B)' in stages[1]['implemented']
    assert 'C=exp(Lambda*realPartSup)' in stages[1]['implemented']
    assert 'axisPressure' not in stages[1]['remaining']
    assert 'H^2 margin' not in stages[1]['remaining']
    assert 'remainderBound/remainderLip' in stages[1]['remaining']
    assert 'realPartSup' in stages[1]['remaining']
    assert 'coefficient-space fixed-point' in stages[1]['remaining']

    assert 'Eq. (5.5)' in stages[2]['implemented']
    assert 'Eq. (5.6)' in stages[2]['implemented']
    assert 'Eq. (5.7)' in stages[2]['implemented']
    assert 'Picard' in stages[2]['implemented']
    assert 'Lemma 5.2' in stages[2]['implemented']
    assert 'five-moment repair' in stages[2]['implemented']
    assert 'Eq. (5.15) forward reconstruction' in stages[2]['implemented']
    assert 'SlowBorelBase/DiagonalScale' in stages[2]['implemented']
    assert 'C[j,m]' in stages[2]['implemented']
    assert 'A0/A1/f_n' in stages[2]['remaining']
    assert 'Eq. (5.6) source' in stages[2]['remaining']
    assert 'true eta-dependent repaired coefficient hierarchy' in stages[2]['remaining']
    assert 'true uniform compactness bounds C[j,m]' in stages[2]['remaining']
    assert 'connect the repaired U_n/E_n through the paper\'s Eq. (5.15)' not in stages[2]['remaining']

    assert 'Eq. (6.8)' in stages[3]['implemented']
    assert 'TangentialBaseJetProvider' in stages[3]['implemented']
    assert 'squared partitions' in stages[3]['implemented']
    assert '2250-color' in stages[3]['implemented']
    assert 'common conservative r0' in stages[3]['implemented']
    assert 'physical slow-support adjacency bridge' in stages[3]['implemented']
    assert 'common physical-point cross-band certificate' in stages[3]['implemented']
    assert 'PhaseEstimates' in stages[3]['implemented']
    assert 'rounded_normal_estimates' in stages[3]['implemented']
    assert 'SlotColoring.Adj/discrete interaction' not in stages[3]['remaining']
    assert 'paper-exact Proposition 5.5 background' in stages[3]['remaining']
    assert 'uniform LocalBaseBounds' in stages[3]['remaining']

    assert 'analytic cutoff gradient' in stages[7]['implemented']
    assert 'time-switch' in stages[7]['implemented']
    assert 'TimeLocalization residual identity' in stages[7]['implemented']
    assert 'closed-past zeroBefore/pastVelocity/pastPressure' in stages[7]['implemented']
    assert 'Taylor-Borel' in stages[7]['implemented']
    assert 'doublingEnvelope' in stages[7]['implemented']
    assert 'full-spacetime endpoint-jet adapter' in stages[7]['implemented']
    assert 'CandidateFromLimits' in stages[7]['implemented']
    assert 'SpatialBorelExtension boundSum/localScale' in stages[7]['implemented']
    assert '2^-j derivative-tail certificates' in stages[7]['implemented']
    assert "closed-past localized NS residual's full spacetime derivative family" in stages[7]['remaining']
    assert 'analytic compact-template derivative bounds' in stages[7]['remaining']
    assert 'smooth global force extension' in stages[7]['remaining']

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
