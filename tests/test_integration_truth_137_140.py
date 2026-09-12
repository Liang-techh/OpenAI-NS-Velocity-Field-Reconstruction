import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_137_140.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_137_140_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "1b7f67397bf6561e7c1f97bdad92099775b33c3e"
    assert ledger["integrated_prs"] == [137, 138, 139, 140]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {137, 138, 139, 140}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["remaining_boundary"].strip()

    assert "complete coefficientOperators record" in rows[137]["remaining_boundary"]
    assert "naturalRemainder(x0)" in rows[137]["remaining_boundary"]
    assert "Section5LowerHistoryJetHierarchy does not yet own this layer" in rows[138]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in rows[138]["remaining_boundary"]
    assert "No paper-exact Proposition 5.5 base-field provider" in rows[139]["remaining_boundary"]
    assert "actual Eq. (9.21) residual sequence" in rows[140]["remaining_boundary"]
    assert "smooth compact forcing" in rows[140]["remaining_boundary"]


def test_source_addenda_canonical_manifest_and_runtime_remain_fail_closed() -> None:
    addenda = (
        "provenance_manifest_addendum_axis_coefficient_inverse_dot_product.json",
        "provenance_manifest_addendum_138.json",
        "provenance_manifest_addendum_phase_large_band_base_source.json",
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


def test_new_source_bindings_do_not_promote_upstream_inputs_to_constructed_fields() -> None:
    stage3 = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_phase_large_band_base_source.json"
    )["layer"]
    stage7 = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_section9_endpoint_residual_source.json"
    )["layer"]

    assert stage3["status"] == "formal-structure"
    assert "No paper-exact Proposition 5.5 provider" in stage3["remaining_boundary"]
    assert "oscillatory wave" in stage3["remaining_boundary"]

    assert stage7["status"] == "formal-structure"
    assert "actual Eq. (9.21) residual sequence" in stage7["remaining_boundary"]
    assert "machine-derived all-order uniform endpoint majorants" in stage7["remaining_boundary"]
    assert "smooth compact forcing" in stage7["remaining_boundary"]
