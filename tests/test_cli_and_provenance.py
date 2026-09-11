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
    assert any('finalAngular' in item and 'clockWeight' in item
               for item in leading['implemented_components'])
    assert any('actual SchedulePressure.axisPressure' in item
               for item in leading['implemented_components'])
    assert any('uniform positive H^2 margin' in item
               for item in leading['missing_for_paper_exact'])
    assert any('coefficient-space fixed-point fields' in item
               for item in leading['missing_for_paper_exact'])

    background=layers['stage-2-all-order-background']
    assert background['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/background_moment_repair.py' in background['artifacts']
    assert any('five-moment repair' in item for item in background['implemented_components'])
    assert any('Eq. (5.15)' in item for item in background['missing_for_paper_exact'])

    oscillatory=layers['stage-3-to-6-oscillatory-corrections']
    assert oscillatory['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/slow_partition.py' in oscillatory['artifacts']
    assert any('squared partitions' in item for item in oscillatory['implemented_components'])
    assert any('Lemma 6.1' in item for item in oscillatory['missing_for_paper_exact'])

    final=layers['stage-7-final-localization']
    assert final['status']=='formal-structure'
    assert 'src/openai_ns_reconstruction/endpoint_jets.py' in final['artifacts']
    assert any('full-spacetime' in item and 'timeVector' in item
               for item in final['implemented_components'])
    assert any('actual localized NS residual' in item
               for item in final['missing_for_paper_exact'])

    assert manifest['full_reconstruction'] is False
