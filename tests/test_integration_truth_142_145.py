import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_142_145.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_142_145_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "c97eeed84c8e03cdb2973e06258e1056afc12f08"
    assert ledger["integrated_prs"] == [142, 143, 144, 145]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {142, 143, 144, 145}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["remaining_boundary"].strip()

    assert "complete coefficientOperators record" in rows[142]["remaining_boundary"]
    assert "naturalRemainder(x0)" in rows[142]["remaining_boundary"]
    assert "hierarchy-owned partial_eta(Omega/X)" in rows[143]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in rows[143]["remaining_boundary"]
    assert "No actual asymptotic partition/carrier/slot" in rows[144]["remaining_boundary"]
    assert "paper-exact Proposition 5.5 base provider" in rows[144]["remaining_boundary"]
    assert "actual Eq. (9.21) residual sequence" in rows[145]["remaining_boundary"]
    assert "smooth compact forcing" in rows[145]["remaining_boundary"]


def test_source_addenda_canonical_manifest_and_runtime_remain_fail_closed() -> None:
    addenda = (
        "provenance_manifest_addendum_axis_coefficient_inverse_mixed.json",
        "provenance_manifest_addendum_143.json",
        "provenance_manifest_addendum_phase_large_band_family_inputs.json",
        "provenance_manifest_addendum_section9_endpoint_residual_source.json",
    )
    for name in addenda:
        data = _load(ROOT / "references" / name)
        assert data["full_reconstruction"] is False
        assert data["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[3]["status"] == "formal-structure"
    assert stages[7]["status"] == "formal-structure"


def test_latest_bindings_preserve_upstream_dependencies_and_evidence_semantics() -> None:
    stage1 = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_inverse_mixed.json"
    )["layer"]
    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_143.json")["layers"][0]
    stage3 = _load(
        ROOT / "references" / "provenance_manifest_addendum_phase_large_band_family_inputs.json"
    )["layer"]
    stage7 = _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_endpoint_residual_source.json"
    )["layer"]

    assert stage1["status"] == "formal-structure"
    assert "complete coefficientOperators record is not yet assembled" in stage1["remaining_boundary"]
    assert "naturalRemainder(x0)" in stage1["remaining_boundary"]

    assert stage2["status"] == "formal-structure"
    assert "order-zero leading data remain upstream" in stage2["remaining_boundary"]
    assert "genuine k=1 Picard step" in stage2["remaining_boundary"]

    assert stage3["status"] == "formal-structure"
    assert "No actual asymptotic partition/carrier/slot" in stage3["remaining_boundary"]
    assert "oscillatory wave" in stage3["remaining_boundary"]

    assert stage7["status"] == "formal-structure"
    assert "ordered-majorant-evidence-kind" in _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_endpoint_residual_source.json"
    )["scope_note"]
    assert "actual Eq. (9.21) residual sequence" in stage7["remaining_boundary"]
    assert "machine-derived all-order uniform endpoint majorants" in stage7["remaining_boundary"]
    assert "smooth compact forcing" in stage7["remaining_boundary"]
