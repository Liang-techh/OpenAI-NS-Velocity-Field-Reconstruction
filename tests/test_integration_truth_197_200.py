import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_197_200.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_197_200_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "72e0d9109823db445d426e96311abeb21e6140d2"
    assert ledger["integrated_prs"] == [197, 198, 199, 200]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {197, 198, 199, 200}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "complete typed mixed-scale naturalRemainder(x0) pair" in rows[197]["capability"]
    assert "1/(2*Lambda)" in rows[197]["remaining_boundary"]
    assert "partial_eta^2(G f_n)" in rows[198]["equation_or_source"]
    assert "caller-provided" in rows[198]["remaining_boundary"]
    assert "one shared Prepared identity" in rows[199]["capability"]
    assert "CellIndex" in rows[199]["remaining_boundary"]
    assert "3/8, 3/4" in rows[200]["capability"]
    assert "smooth extension through t=1" in rows[200]["remaining_boundary"]


def test_landed_provenance_boundaries_remain_fail_closed() -> None:
    stage1 = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_natural_remainder.json"
    )
    assert stage1["full_reconstruction"] is False
    assert stage1["paper_exact_velocity_available"] is False
    assert stage1["layer"]["status"] == "formal-structure"
    assert "complete mixed-scale naturalRemainder(x0) pair" in stage1["layer"]["capability"]
    assert "1/(2*Lambda)" in stage1["layer"]["remaining_boundary"]
    assert "genuine Picard x1" in stage1["layer"]["remaining_boundary"]

    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_198.json")
    assert stage2["full_reconstruction"] is False
    assert stage2["paper_exact_velocity_available"] is False
    assert stage2["layers"][0]["status"] == "formal-structure"
    assert "explicit caller data" in stage2["layers"][0]["remaining_boundary"]
    assert "Picard iteration/convergence" in stage2["layers"][0]["remaining_boundary"]

    stage3 = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_CANONICAL_SCOPE_COVERAGE_PROVENANCE.json"
    )
    assert stage3["full_reconstruction"] is False
    assert stage3["paper_exact_velocity_available"] is False
    truth = stage3["layer"]["truth_boundary"]
    assert truth["actual_canonical_prepared_application_machine_verified"] is False
    assert truth["actual_cell_index_enumeration_machine_verified"] is False
    assert truth["actual_physical_base_values_materialized"] is False
    assert truth["uniform_eq_7_9_to_7_11_verified"] is False
    assert truth["paper_exact_velocity_available"] is False

    section10 = (
        ROOT / "references" / "SECTION10_PAPER_TIME_SWITCH_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "does not construct the missing one-sided Section 9 field extension through `t = 1`" in section10
    assert "does not promote paper-exact velocity or forcing" in section10
    assert "section9_field_smooth_extension_through_t1_constructed`" in section10
    assert "paper_exact_velocity_available`" in section10

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["full_paper_reconstruction_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_are_explicit() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any(
        "naturalRemainder(x0)" in item
        and "1/(2*Lambda)" in item
        and "fixed-point" in item
        for item in forbidden
    )
    assert any(
        "partial_eta^2 f_n" in item
        and "Picard convergence" in item
        and "Proposition 5.3" in item
        for item in forbidden
    )
    assert any(
        "Prepared" in item
        and "CellIndex" in item
        and "Eqs. (7.9)-(7.11)" in item
        for item in forbidden
    )
    assert any(
        "timeSwitch" in item
        and "t=1" in item
        and "final forcing" in item
        for item in forbidden
    )
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)
