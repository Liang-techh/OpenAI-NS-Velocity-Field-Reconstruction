import json

from openai_ns_reconstruction.cli import construction_audit, main
from openai_ns_reconstruction.status import construction_status


def test_audit_derives_truth_from_status_without_promotion():
    status = construction_status()
    audit = construction_audit()

    assert audit["kind"] == "construction-audit"
    assert audit["paper_exact_velocity_available"] is False
    assert audit["completion_gate_passed"] is False
    assert audit["paper_exact_velocity_available"] == status["paper_exact_velocity_available"]
    assert [stage["status"] for stage in audit["stages"]] == [
        stage["status"] for stage in status["stages"]
    ]
    assert [stage["remaining"] for stage in audit["stages"]] == [
        stage["remaining"] for stage in status["stages"]
    ]


def test_audit_cli_exists_and_paper_exact_gate_fails_closed(capsys):
    rc = main(["audit", "--require-paper-exact"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert rc == 2
    assert payload["kind"] == "construction-audit"
    assert payload["paper_exact_velocity_available"] is False
    assert payload["completion_gate_passed"] is False
    assert "Incomplete paper construction" in captured.err


def test_demo_paper_exact_gate_refuses_before_writing(tmp_path, capsys):
    output = tmp_path / "must-not-exist"
    rc = main(["demo", "--output", str(output), "--require-paper-exact"])
    captured = capsys.readouterr()

    assert rc == 2
    assert captured.out == ""
    assert "Incomplete paper construction" in captured.err
    assert not output.exists()
