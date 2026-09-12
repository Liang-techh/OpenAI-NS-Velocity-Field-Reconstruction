import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_193_195.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_193_195_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "f02883fb0a9a298857ad3190e5783063f260f38a"
    assert ledger["integrated_prs"] == [193, 194, 195]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {193, 194, 195}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        if row["provenance"] is not None:
            assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "mixed-scale angular branch" in rows[193]["capability"]
    assert "complete naturalRemainder(x0)" in rows[193]["remaining_boundary"]
    assert "canonicalPrepared/canonicalPhases" in rows[194]["capability"]
    assert "CellIndex" in rows[194]["remaining_boundary"]
    assert "Eq. (4.1)" in rows[195]["capability"]
    assert "cutoff derivative table is still provider input" in rows[195]["remaining_boundary"]


def test_landed_source_truth_boundaries_remain_fail_closed() -> None:
    angular = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_angular_remainder.json"
    )
    assert angular["full_reconstruction"] is False
    assert angular["paper_exact_velocity_available"] is False
    assert angular["layer"]["status"] == "formal-structure"
    assert angular["layer"]["actual_data"]["fixed_point_state"] == (
        "x0 = actual SchedulePressure referencePair"
    )
    assert "angular_formula" in angular["layer"]["actual_data"]
    assert "complete naturalRemainder(x0)" in angular["layer"]["remaining_boundary"]
    assert "Picard x1" in angular["layer"]["remaining_boundary"]

    canonical_base = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_CANONICAL_PHYSICAL_BASE_PROVENANCE.json"
    )
    assert canonical_base["full_reconstruction"] is False
    assert canonical_base["paper_exact_velocity_available"] is False
    assert canonical_base["layer"]["status"] == "formal-structure"
    truth = canonical_base["layer"]["truth_boundary"]
    assert truth["actual_canonical_prepared_application_machine_verified"] is False
    assert truth["actual_cell_index_enumeration_machine_verified"] is False
    assert truth["actual_physical_base_values_materialized"] is False
    assert truth["actual_base_fields_verified"] is False
    assert truth["uniform_eq_7_9_to_7_11_verified"] is False
    assert truth["paper_exact_velocity_available"] is False

    eq41_source = (
        ROOT / "src" / "openai_ns_reconstruction" / "section9_eq41_similarity_jet.py"
    ).read_text(encoding="utf-8")
    assert "No finite differences, interpolation, fitting, or" in eq41_source
    assert "one-dimensional cutoff derivative table with the manuscript's fixed cutoff" in eq41_source
    assert "does not construct the infinite Eq. (9.21) sum or its t=1 extension" in eq41_source
    assert "paper_exact_velocity_available: bool = False" in eq41_source

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_are_explicit() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any(
        "angular and axial" in item
        and "complete naturalRemainder(x0)" in item
        and "fixed point" in item
        for item in forbidden
    )
    assert any(
        "canonicalPrepared/canonicalPhases" in item
        and "CellIndex" in item
        and "FinalSlowBase.velocity" in item
        for item in forbidden
    )
    assert any(
        "Eq. (4.1)" in item
        and "provider-supplied cutoff derivatives" in item
        and "infinite Eq. (9.21)" in item
        for item in forbidden
    )
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)
