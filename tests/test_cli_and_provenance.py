import hashlib
import json
import pytest
from openai_ns_reconstruction.__main__ import main
from openai_ns_reconstruction.provenance import (
    status_report,require_complete_reconstruction,IncompleteReconstructionError,
)
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
    assert not status_report()['sources']['lean_build_verified']
    assert status_report()['stages'][0]['complete']


def test_exact_audit_returns_nonzero_without_fake_success(tmp_path,capsys):
    path=tmp_path/'audit.json'
    assert main(['audit','--output',str(path),'--require-paper-exact'])==2
    assert json.loads(path.read_text())['full_reconstruction'] is False


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
    import json
    from pathlib import Path
    from openai_ns_reconstruction.provenance import status_report
    manifest = json.loads((Path(__file__).resolve().parents[1] /
                           "references/provenance_manifest.json").read_text())
    report = status_report()
    assert manifest["sources"]["official_lean"]["commit"] == report["sources"]["upstream_commit"]
    assert manifest["sources"]["paper"]["url"] == report["sources"]["paper_url"]
    assert manifest["sources"]["paper"]["sha256"] == report["sources"]["paper_sha256"]
    assert manifest["sources"]["official_lean"]["build_verified"] is False
    assert manifest["full_reconstruction"] == report["full_reconstruction"] is False
