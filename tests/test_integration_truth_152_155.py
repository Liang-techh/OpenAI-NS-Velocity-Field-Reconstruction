import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_152_155.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_152_155_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "1766ceb9d33feda4baf255f05679ff3c280435ba"
    assert ledger["integrated_prs"] == [152, 153, 154, 155]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {152, 153, 154, 155}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["remaining_boundary"].strip()

    assert "naturalRemainder(x0)" in rows[152]["remaining_boundary"]
    assert "first genuine Picard iterate" in rows[152]["remaining_boundary"]
    assert "hierarchy-owned partial_eta actualLowerSource" in rows[153]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in rows[153]["remaining_boundary"]
    assert "paper-exact Proposition 5.5 base/background provider" in rows[154]["remaining_boundary"]
    assert "uniform Eq. (7.9)-(7.11) verification" in rows[154]["remaining_boundary"]
    assert "actual Eq. (9.21) residual sequence" in rows[155]["remaining_boundary"]
    assert "smooth compact forcing" in rows[155]["remaining_boundary"]


def test_source_provenance_canonical_manifest_and_runtime_remain_fail_closed() -> None:
    stage1 = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_data.json"
    )
    assert stage1["full_reconstruction"] is False
    assert stage1["paper_exact_velocity_available"] is False
    assert stage1["layer"]["status"] == "formal-structure"

    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_153.json")
    assert stage2["full_reconstruction"] is False
    assert stage2["paper_exact_velocity_available"] is False
    assert stage2["layers"][0]["status"] == "formal-structure"

    stage3 = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_COORDINATE_ERRORS_PROVENANCE.json"
    )
    assert stage3["status"] == "formal-structure"
    assert stage3["truth_flags"]["actual_coordinate_errors_materialized"] is False
    assert stage3["truth_flags"]["uniform_eq_7_9_to_7_11_verified"] is False
    assert stage3["truth_flags"]["paper_exact_velocity_available"] is False

    stage7_text = (
        ROOT / "references" / "SECTION9_ENDPOINT_MAJORANT_SOURCE_BINDING_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "source_majorants_derived_from_actual_residual_verified" in stage7_text
    assert "actual_section9_sequence_verified" in stage7_text
    assert "endpoint-limit construction/uniqueness" in stage7_text
    assert "paper_exact_velocity_available" in stage7_text

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


def test_latest_capabilities_do_not_erase_upstream_or_endpoint_boundaries() -> None:
    stage1 = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_data.json"
    )["layer"]
    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_153.json")["layers"][0]
    stage3 = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_COORDINATE_ERRORS_PROVENANCE.json"
    )
    stage7_text = (
        ROOT / "references" / "SECTION9_ENDPOINT_MAJORANT_SOURCE_BINDING_PROVENANCE.md"
    ).read_text(encoding="utf-8")

    assert "naturalRemainder composition has not yet been executed" in stage1["remaining_boundary"]
    assert "first genuine Picard iterate" in stage1["remaining_boundary"]

    assert "hierarchy-owned partial_eta actualLowerSource" in stage2["remaining_boundary"]
    assert "genuine k=1 Eq. (5.7) Picard application" in stage2["remaining_boundary"]

    remaining_stage3 = " ".join(stage3["remaining_for_paper_exact"])
    assert "genuine Proposition 5.5 base/background provider" in remaining_stage3
    assert "curl-realized oscillatory waves" in remaining_stage3

    assert "does **not** derive the uniform Eq. (9.18) constants" in stage7_text
    assert "actual Eq. (9.21) residual construction" in stage7_text
    assert "infinite Borel right-jet convergence" in stage7_text
    assert "smooth compact forcing" in stage7_text
    assert "bounded-energy and blow-up closure" in stage7_text
