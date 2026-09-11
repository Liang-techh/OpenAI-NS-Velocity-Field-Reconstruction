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
    assert 'AxisContraction-style bound/Lipschitz propagation' in stages[1]['implemented']
    assert 'radiusLoss(1/2)=12' in stages[1]['implemented']
    assert 'normalized amplitude M=12' in stages[1]['implemented']
    assert 'K<=30720B' in stages[1]['implemented']
    assert 'factorial-series enclosure' in stages[1]['implemented']
    assert 'actual-schedule analytic-neighborhood certificate' in stages[1]['implemented']
    assert 'positive common rho' in stages[1]['implemented']
    assert 'eleven-field common complex bound B' in stages[1]['implemented']
    assert 'axisPhase real-part supremum' in stages[1]['implemented']
    assert '96-digit Decimal' in stages[1]['implemented']
    assert 'sys.float_info.max' in stages[1]['implemented']
    assert 'landed wide Lambda/C selector' in stages[1]['implemented']
    assert 'symbolic exponent' in stages[1]['implemented']
    assert 'without binary64 down-conversion' in stages[1]['implemented']
    assert 'axis_fixed_point_picard' in stages[1]['implemented']
    assert 's*remainderBound<=1' in stages[1]['implemented']
    assert 's*remainderLip<=1/2' in stages[1]['implemented']
    assert 'axis_reference_pair' in stages[1]['implemented']
    assert 'zeroth parameter-jet radial coefficient functions' in stages[1]['implemented']
    assert 'axisPressure' not in stages[1]['remaining']
    assert 'H^2 margin' not in stages[1]['remaining']
    assert 'actual admissible common analytic radius' not in stages[1]['remaining']
    assert 'common complex-field sup bound' not in stages[1]['remaining']
    assert 'AxisCoefficientSpace' in stages[1]['remaining']
    assert 'coefficientOperators/naturalRemainder' in stages[1]['remaining']
    assert 'materialize the coefficient-space fixed-point phi/u' in stages[1]['remaining']
    assert 'Derive the corresponding average/pressure fields' in stages[1]['remaining']
    assert 'former binary64 Lambda/C representation blocker is closed' in stages[1]['remaining']

    assert 'Eq. (5.5)' in stages[2]['implemented']
    assert 'Eq. (5.6)' in stages[2]['implemented']
    assert 'Eq. (5.7)' in stages[2]['implemented']
    assert 'Picard' in stages[2]['implemented']
    assert 'PositiveAxisSystem adapter' in stages[2]['implemented']
    assert 'BaseJet(phi,axial,beta)' in stages[2]['implemented']
    assert 'SourceJet(angular,axial,pressureProduct,omegaQuotient)' in stages[2]['implemented']
    assert 'lower-history bridge' in stages[2]['implemented']
    assert 'strict lower-order ProfileSecondJet history' in stages[2]['implemented']
    assert 'Eq. (5.8)' in stages[2]['implemented']
    assert 'p_k=ceil(k/2)' in stages[2]['implemented']
    assert 'complete-tail truncation certificate' in stages[2]['implemented']
    assert 'Lemma 5.2' in stages[2]['implemented']
    assert 'five-moment repair' in stages[2]['implemented']
    assert 'function-level compact-repair adapter' in stages[2]['implemented']
    assert 'paper_exact=False' in stages[2]['implemented']
    assert 'Eq. (5.15) forward reconstruction' in stages[2]['implemented']
    assert 'SlowBorelBase/DiagonalScale' in stages[2]['implemented']
    assert 'C[j,m]' in stages[2]['implemented']
    assert 'arbitrary-precision power-of-two Python integer witness' in stages[2]['implemented']
    assert 'reciprocal_support_log_edge' in stages[2]['implemented']
    assert 'finite-prefix cutoff support certificate' in stages[2]['implemented']
    assert 'unique possible transition order' in stages[2]['implemented']
    assert 'stable truncation' in stages[2]['implemented']
    assert '2^-J q^(h(J+1)-m)' in stages[2]['implemented']
    assert 'h(J+1)+b-2M' in stages[2]['implemented']
    assert 'minimal J' in stages[2]['implemented']
    assert 'A0/A1/f_n' in stages[2]['remaining']
    assert 'actual C_n' in stages[2]['remaining']
    assert 'Eq. (5.6) source' in stages[2]['remaining']
    assert 'true eta-dependent repaired hierarchy' in stages[2]['remaining']
    assert 'true uniform C[j,m] bounds' in stages[2]['remaining']
    assert 'theorem-level local finiteness' in stages[2]['remaining']
    assert 'all-jets-flat residual decay' in stages[2]['remaining']
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
    assert 'UniformLocalBase bridge' in stages[3]['implemented']
    assert 'd<=S^-3' in stages[3]['implemented']
    assert 'M(S^-3+epsilon^2)' in stages[3]['implemented']
    assert 'BasePhaseGeometry' in stages[3]['implemented']
    assert 'phaseConstant(M)=normalConstant(frequencyBound(M))' in stages[3]['implemented']
    assert '16 G^2(1+G)E' in stages[3]['implemented']
    assert '4 M(2A+5)delta' in stages[3]['implemented']
    assert 'SlotColoring.Adj/discrete interaction' not in stages[3]['remaining']
    assert 'paper-exact Proposition 5.5 background' in stages[3]['remaining']
    assert 'actual LocalBaseBounds/C1-C2 hypotheses' in stages[3]['remaining']
    assert 'normal-closeness' in stages[3]['remaining']
    assert 'landed uniform normal and frame/damping envelopes' in stages[3]['remaining']

    assert stages[4]['status']=='formal-structure'
    assert 'signed two-slot reference covariance algebra' in stages[4]['implemented']
    assert 'strict positive-cone condition |a t|<b m' in stages[4]['implemented']
    assert 'epsilon*mask^2 covariance scaling' in stages[4]['implemented']
    assert 'coefficient-level curl-realization algebra' in stages[4]['implemented']
    assert 'cylindrical_curl_jet' in stages[4]['implemented']
    assert 'B_theta/r' in stages[4]['implemented']
    assert 'primary amplitude ODE algebra' in stages[4]['implemented']
    assert 'exact 2x2 modal operator' in stages[4]['implemented']
    assert 'x=p+q, y=h(p-q)' in stages[4]['implemented']
    assert 'finite-interval PrimaryODE adapter' in stages[4]['implemented']
    assert 'Volterra integral defect' in stages[4]['implemented']
    assert 'closed-form constant-coefficient oracle' in stages[4]['implemented']
    assert 'actual pulse-integrated covariance columns' in stages[4]['remaining']
    assert 'determinant/inverse stability' in stages[4]['remaining']
    assert 'finite-interval PrimaryODE Volterra solution adapter' in stages[4]['remaining']
    assert 'rigorous error certificate' in stages[4]['remaining']
    assert 'cylindrical derivative jet accepted by cylindrical_curl_jet' in stages[4]['remaining']
    assert 'genuine supported divergence-free oscillatory wave' in stages[4]['remaining']

    assert 'analytic cutoff gradient' in stages[7]['implemented']
    assert 'time-switch' in stages[7]['implemented']
    assert 'TimeLocalization residual identity' in stages[7]['implemented']
    assert 'closed-past zeroBefore/pastVelocity/pastPressure' in stages[7]['implemented']
    assert 'Taylor-Borel' in stages[7]['implemented']
    assert 'doublingEnvelope' in stages[7]['implemented']
    assert 'close_left_open_past' in stages[7]['implemented']
    assert 'rejects t>T' in stages[7]['implemented']
    assert 'full-spacetime endpoint-jet adapter' in stages[7]['implemented']
    assert 'CandidateFromLimits' in stages[7]['implemented']
    assert 'SpatialBorelExtension' in stages[7]['implemented']
    assert '2^-j derivative-tail certificates' in stages[7]['implemented']
    assert 'CandidateFromLimits trace/glue bridge' in stages[7]['implemented']
    assert 'degree-zero endpoint tensor' in stages[7]['implemented']
    assert 'volume pi/32' in stages[7]['implemented']
    assert 'E(t)<=pi M^2/64' in stages[7]['implemented']
    assert 'late-origin preservation certificate' in stages[7]['implemented']
    assert '3/4<=t<1' in stages[7]['implemented']
    assert 'c=1 and grad(c)=0 at the origin' in stages[7]['implemented']
    assert 'endpoint-localization transfer' in stages[7]['implemented']
    assert 'requires T=1' in stages[7]['implemented']
    assert 'timeSwitch\'=0' in stages[7]['implemented']
    assert 'section10_endpoint_ladder' in stages[7]['implemented']
    assert 'degrees 0..N' in stages[7]['implemented']
    assert 'pointwise endpoint-force support certificate' in stages[7]['implemented']
    assert 'requiring exact zeros with no tolerance' in stages[7]['implemented']
    assert "closed-past localized NS residual's full spacetime derivative family" in stages[7]['remaining']
    assert 'official t>=3/4 plateau' in stages[7]['remaining']
    assert 'close the infinite all-order family' in stages[7]['remaining']
    assert 'degree-zero endpoint value' in stages[7]['remaining']
    assert 'every genuine endpoint normal coefficient vanish outside' in stages[7]['remaining']
    assert 'analytic compact-template derivative bounds' in stages[7]['remaining']
    assert 'smooth through t=1' in stages[7]['remaining']
    assert 'uniform bounded kinetic-energy' in stages[7]['remaining']
    assert '||curl A(t,0)||->infinity premise' in stages[7]['remaining']

    limitations=' '.join(status['limitations'])
    assert 'traced-residual adapters' in limitations
    assert 'caller-supplied endpoint jets' in limitations
    assert 'close_left_open_past' in limitations
    assert 'official-plateau transfer' in limitations
    assert 'finite section10_endpoint_ladder' in limitations
    assert 'endpoint point-support certificate' in limitations
    assert 'finite collection of point queries' in limitations
    assert 'fixed-time implication' in limitations
    assert 'uniform bounded-energy estimate' in limitations
    assert 'late-origin adapter preserves an upstream-certified origin curl value' in limitations
    assert 'primary-amplitude ODE' in limitations
    assert 'generic finite-interval PrimaryODE/Volterra execution adapter' in limitations
    assert 'paper-exact wave are still missing' in limitations
    assert 'wide Lambda/C selector' in limitations
    assert 'landed Picard gate' in limitations
    assert 'referencePair radial coefficient functions' in limitations
    assert 'do not materialize the coefficient-space fixed point' in limitations
    assert stages[1]['status']==stages[2]['status']==stages[3]['status']==stages[4]['status']==stages[7]['status']=='formal-structure'
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
