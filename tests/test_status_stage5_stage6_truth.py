import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_runtime_status_exposes_landed_stage5_and_stage6_formal_structure():
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False

    assert stages[5]["status"] == "formal-structure"
    assert "five-row mean-rank update" in stages[5]["implemented"]
    assert "synthetic/caller-supplied" in stages[5]["implemented"]
    assert "actual mean debt" in stages[5]["remaining"]
    assert "paper-exact mean correction" in stages[5]["remaining"]

    assert stages[6]["status"] == "formal-structure"
    assert "Eqs. (9.17)-(9.18)" in stages[6]["implemented"]
    assert "finite Eq. (9.21) prefix" in stages[6]["implemented"]
    assert "paper_exact_velocity_available false" in stages[6]["implemented"]
    assert "genuine Section 9 correction fields" in stages[6]["remaining"]
    assert "infinite Eq. (9.21) sums" in stages[6]["remaining"]
    assert "Proposition 9.9" in stages[6]["remaining"]


def test_stage5_stage6_status_provenance_stays_fail_closed():
    canonical = json.loads((ROOT / "references/provenance_manifest.json").read_text())
    addendum = json.loads(
        (ROOT / "references/provenance_manifest_addendum_status_stage5_6.json").read_text()
    )

    assert canonical["full_reconstruction"] is False
    layer = next(
        layer
        for layer in canonical["layers"]
        if layer["id"] == "stage-3-to-6-oscillatory-corrections"
    )
    assert layer["status"] == "formal-structure"

    assert addendum["canonical_truth_required"]["full_reconstruction"] is False
    assert addendum["canonical_truth_required"]["paper_exact_velocity_available"] is False
    assert addendum["status_surface_change"]["promotion_ceiling"] == "formal-structure"
    assert addendum["status_surface_change"]["paper_exact_gate_changed"] is False

    staged = {entry["stage"]: entry for entry in addendum["stages"]}
    assert staged[5]["status"] == staged[6]["status"] == "formal-structure"
    assert any("actual mean debt" in item for item in staged[5]["not_verified"])
    assert any("infinite Eq. (9.21)" in item for item in staged[6]["not_verified"])
    assert any("endpoint smooth forcing" in item for item in staged[6]["not_verified"])

    for entry in addendum["stages"]:
        for artifact in entry["artifacts"]:
            assert (ROOT / artifact).exists(), artifact
