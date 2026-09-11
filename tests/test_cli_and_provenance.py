import hashlib
import json
import pytest
from openai_ns_reconstruction.__main__ import main
from openai_ns_reconstruction.provenance import (
    status_report,require_complete_reconstruction,IncompleteReconstructionError,
)
from openai_ns_reconstruction.status import construction_status
from openai_ns_reconstruction.demo import toy_gaussian_profile


def test_incomplete_reconstruction_fails_closed():
    report=status_report()
    assert not report['full_reconstruction']
    assert report['blockers']
    assert not report['sources']['lean_build_verified']
    with pytest.raises(IncompleteReconstructionError):
        require_complete_reconstruction()


def test_status_report_has_no_shared_mutable_state():
    report=status_report()
    report['sources']['lean_build_verified']=True
    report['stages'][0]['complete']=False
    report['stages'][1]['remaining']='fake'
    assert not status_report()['sources']['lean_build_verified']
    assert status_report()['stages'][0]['complete']
    assert status_report()['stages'][1]['remaining']!='fake'


def test_audit_and_status_share_one_truth_surface():
    audit=status_report()
    runtime=construction_status()
    assert audit['full_reconstruction'] is runtime['paper_exact_velocity_available'] is False
    audit_stages={stage['stage']:stage for stage in audit['stages']}
    runtime_stages={stage['id']:stage for stage in runtime['stages']}
    assert audit_stages.keys()==runtime_stages.keys()
    for stage_id, stage in runtime_stages.items():
        assert audit_stages[stage_id]['name']==stage['name']
        assert audit_stages[stage_id]['status']==stage['status']
        if 'implemented' in stage:
            assert audit_stages[stage_id]['implemented']==stage['implemented']
        if 'remaining' in stage:
            assert audit_stages[stage_id]['remaining']==stage['remaining']


def test_exact_audit_returns_nonzero_without_fake_success(tmp_path,capsys):
    path=tmp_path/'audit.json'
    assert main(['audit','--output',str(path),'--require-paper-exact'])==2
    report=json.loads(path.read_text())
    assert report['full_reconstruction'] is False
    assert report['paper_exact_velocity_available'] is False


def test_exact_demo_refuses_before_writing_toy_artifacts(tmp_path,capsys):
    path=tmp_path/'not-created'
    assert main(['demo','--output',str(path),'--require-paper-exact'])==2
    assert not path.exists()


def test_demo_runs_offline_and_hashes_artifacts(tmp_path,capsys):
    assert main(['demo','--output',str(tmp_path)])==0
    report=json.loads((tmp_path/'report.json').read_text())
    assert report['all_diagnostic_checks_passed']
    assert report['profile_status']=='toy'
    assert report['full_reconstruction'] is False
    for filename,digest in report['artifacts'].items():
        assert hashlib.sha256((tmp_path/filename).read_bytes()).hexdigest()==digest
    assert len(report['artifacts'])==3
    assert toy_gaussian_profile().paper_exact is False


def test_runtime_source_pins_match_repository_manifest():
    from pathlib import Path
    manifest = json.loads((Path(__file__).resolve().parents[1] /
                           "references/provenance_manifest.json").read_text())
    report = status_report()
    assert manifest["sources"]["official_lean"]["commit"] == report["sources"]["upstream_commit"]
    assert manifest["sources"]["paper"]["url"] == report["sources"]["paper_url"]
    assert manifest["sources"]["paper"]["sha256"] == report["sources"]["paper_sha256"]
    assert manifest["sources"]["official_lean"]["build_verified"] is False
    assert manifest["full_reconstruction"] == report["full_reconstruction"] is False


def test_manifest_tracks_landed_formal_structure_without_promoting_completion():
    from pathlib import Path
    manifest = json.loads((Path(__file__).resolve().parents[1] /
                           "references/provenance_manifest.json").read_text())
    layers={layer['id']:layer for layer in manifest['layers']}

    leading=layers['stage-1-leading-profile']
    assert leading['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/outgoing_tail.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/schedule_axis_pressure.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/schedule_axis_margin.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/schedule_analytic_neighborhood.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/natural_scale_selection.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/natural_scale_selection_wide.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/axis_remainder_bounds.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/axis_remainder_wide_bounds.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/axis_analytic_input_bounds.py' in leading['artifacts']
    assert 'src/openai_ns_reconstruction/stage1_scale_chain_wide.py' in leading['artifacts']
    assert any('finalAngular' in item and 'clockWeight' in item
               for item in leading['implemented_components'])
    assert any('actual SchedulePressure.axisPressure' in item
               for item in leading['implemented_components'])
    assert any('low-|Z|' in item and 'sigma=sqrt(m)/20' in item
               for item in leading['implemented_components'])
    assert any('Lambda=max(1+B+L,1+(14000/9)B)' in item and 'C=exp' in item
               for item in leading['implemented_components'])
    assert any('remainderBound/remainderLip propagation' in item
               for item in leading['implemented_components'])
    assert any('radiusLoss(1/2)=12' in item and 'M=12' in item and 'K<=30720B' in item
               for item in leading['implemented_components'])
    assert any('actual-schedule analytic-neighborhood certificate' in item and 'realPartSup' in item
               for item in leading['implemented_components'])
    assert any('96-digit Decimal' in item and 'sys.float_info.max' in item
               for item in leading['implemented_components'])
    assert any('wide theorem-shaped Lambda/C continuation' in item
               and 'SymbolicExponentialThreshold' in item
               and 'without binary64 down-conversion' in item
               for item in leading['implemented_components'])
    assert not any('materialize and certify an actual admissible common analytic' in item
                   for item in leading['missing_for_paper_exact'])
    assert not any('wide Decimal remainderBound/remainderLip' in item and 'through the pinned Lambda/C' in item
                   for item in leading['missing_for_paper_exact'])
    assert any('coefficient-space fixed-point fields phi/u' in item
               for item in leading['missing_for_paper_exact'])
    assert not any('naturalResolvent factorial-series norm bound' in item
                   for item in leading['missing_for_paper_exact'])

    background=layers['stage-2-all-order-background']
    assert background['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/background_picard_bounds.py' in background['artifacts']
    assert 'src/openai_ns_reconstruction/background_moment_repair.py' in background['artifacts']
    assert 'src/openai_ns_reconstruction/background_moment_repair_profile.py' in background['artifacts']
    assert 'src/openai_ns_reconstruction/background_extension.py' in background['artifacts']
    assert 'src/openai_ns_reconstruction/background_cutoff_schedule.py' in background['artifacts']
    assert 'src/openai_ns_reconstruction/background_cutoff_support.py' in background['artifacts']
    assert 'src/openai_ns_reconstruction/background_tail_order.py' in background['artifacts']
    assert any('five-moment repair' in item for item in background['implemented_components'])
    assert any('function-level Lemma 5.2 compact repair' in item and 'paper_exact' in item
               for item in background['implemented_components'])
    assert any('Eq. (5.8)' in item and 'p_k=ceil(k/2)' in item and 'complete-tail' in item
               for item in background['implemented_components'])
    assert any('Eq. (5.15) forward reconstruction' in item
               for item in background['implemented_components'])
    assert any('SlowBorelBase/DiagonalScale' in item and 'C[j,m]' in item
               for item in background['implemented_components'])
    assert any('arbitrary-precision power-of-two integer witness' in item
               and 'reciprocal_support_log_edge' in item
               for item in background['implemented_components'])
    assert any('finite-prefix cutoff support/plateau classifier' in item and 'stable truncation' in item
               for item in background['implemented_components'])
    assert any('2^-J' in item and 'h(J+1)-m' in item and 'h(J+1)+b-2M' in item
               for item in background['implemented_components'])
    assert any('A0/A1/f_n' in item and 'C_n' in item and 'Delta' in item
               for item in background['missing_for_paper_exact'])
    assert any('true eta-dependent repaired coefficient hierarchy' in item
               for item in background['missing_for_paper_exact'])
    assert any('true uniform compactness bounds C[j,m]' in item and 'theorem-level local finiteness' in item
               for item in background['missing_for_paper_exact'])
    assert any('all-jets-flat residual decay' in item for item in background['missing_for_paper_exact'])

    oscillatory=layers['stage-3-to-6-oscillatory-corrections']
    assert oscillatory['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/slow_partition.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/slot_geometry.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/slow_support_adjacency.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/phase_uniform_bounds.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/phase_frame_bounds.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/stress_cone.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/pulse_covariance_budget.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/curl_realization_algebra.py' in oscillatory['artifacts']
    assert 'src/openai_ns_reconstruction/primary_amplitude_ode.py' in oscillatory['artifacts']
    assert any('squared partitions' in item for item in oscillatory['implemented_components'])
    assert any('2250-color' in item and 'common r0' in item
               for item in oscillatory['implemented_components'])
    assert any('physical slow-support to SlotColoring.Adj bridge' in item
               for item in oscillatory['implemented_components'])
    assert any('UniformLocalBase bridge' in item and 'M(S^-3+epsilon^2)' in item
               for item in oscillatory['implemented_components'])
    assert any('phaseConstant(M)=normalConstant(frequencyBound(M))' in item
               and '16 G^2(1+G)E' in item and '4 M(2A+5)delta' in item
               for item in oscillatory['implemented_components'])
    assert any('strict positive-cone test |a t|<b m' in item and 'sqrt(epsilon)*mask' in item
               for item in oscillatory['implemented_components'])
    assert any('Formula (30)' in item and 'inverseCarrier=i/K' in item
               for item in oscillatory['implemented_components'])
    assert any('primary amplitude ODE pointwise algebra' in item and 'modalOperator' in item
               for item in oscillatory['implemented_components'])
    assert any('paper-exact Proposition 5.5 background' in item
               for item in oscillatory['missing_for_paper_exact'])
    assert any('LocalBaseBounds C1/C2 hypotheses' in item and 'normal-closeness' in item
               for item in oscillatory['missing_for_paper_exact'])
    assert any('actual pulse-integrated covariance matrix H' in item and 'Eq. (7.28)' in item
               for item in oscillatory['missing_for_paper_exact'])
    assert any('finite-interval PrimaryODE Volterra' in item and 'LocalizedCurlRealization' in item
               for item in oscillatory['missing_for_paper_exact'])

    final=layers['stage-7-final-localization']
    assert final['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/past_extension.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/endpoint_jets.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/endpoint_scale_schedule.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/endpoint_limit_majorant.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/endpoint_support.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/traced_residual.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/section10_support_energy.py' in final['artifacts']
    assert 'src/openai_ns_reconstruction/section10_origin_preservation.py' in final['artifacts']
    assert any('PastExtension closed-past branch' in item
               for item in final['implemented_components'])
    assert any('close_left_open_past' in item and 'rejects t>T' in item
               for item in final['implemented_components'])
    assert any('full-spacetime' in item and 'timeVector' in item
               for item in final['implemented_components'])
    assert any('boundSum/localScale' in item and '2^-j' in item
               for item in final['implemented_components'])
    assert any('CandidateFromLimits traced-residual/Borel-glue bridge' in item
               for item in final['implemented_components'])
    assert any('pi/32' in item and 'E(t)<=pi M^2/64' in item
               for item in final['implemented_components'])
    assert any('late-origin localization certificate' in item and '3/4<=t<1' in item
               for item in final['implemented_components'])
    assert any('pointwise endpoint-force support implication' in item and 'no tolerance' in item
               for item in final['implemented_components'])
    assert any('closed-past localized NS residual' in item and 'full spacetime derivative family' in item
               for item in final['missing_for_paper_exact'])
    assert any('degree-zero endpoint value used by close_left_open_past' in item
               for item in final['missing_for_paper_exact'])
    assert any('analytic compact-template derivative bounds' in item
               for item in final['missing_for_paper_exact'])
    assert any('global compact-support conclusion' in item
               for item in final['missing_for_paper_exact'])
    assert any('uniform bounded kinetic-energy' in item
               for item in final['missing_for_paper_exact'])
    assert any('SpatialCurl.spatialCurl A' in item and 'late-origin certificate' in item
               for item in final['missing_for_paper_exact'])
    assert any('tautological' in item for item in final['numerical_caveats'])
    assert any('supplied endpoint tensors' in item and 'not evidence' in item
               for item in final['numerical_caveats'])
    assert any('endpoint_borel.close_left_open_past' in item and 'does not prove continuity' in item
               for item in final['numerical_caveats'])
    assert any('endpoint_support.py' in item and 'finitely many point queries' in item
               for item in final['numerical_caveats'])
    assert any('section10_support_energy.py' in item and 'fixed-time implication' in item
               for item in final['numerical_caveats'])
    assert any('section10_origin_preservation.py' in item and 'not evidence' in item
               for item in final['numerical_caveats'])

    runtime = construction_status()
    runtime_stages = {stage['id']: stage for stage in runtime['stages']}
    assert runtime_stages[4]['status'] == 'formal-structure'
    assert 'strict positive-cone condition' in runtime_stages[4]['implemented']
    assert 'primary amplitude ODE algebra' in runtime_stages[4]['implemented']
    assert 'actual pulse-integrated covariance' in runtime_stages[4]['remaining']
    assert 'finite-interval PrimaryODE Volterra solution' in runtime_stages[4]['remaining']
    assert runtime['paper_exact_velocity_available'] is False
    assert manifest['full_reconstruction'] is False
