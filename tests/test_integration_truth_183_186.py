import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_183_186.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_183_186_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "4ff3f85634a32fc706a92dc0caf886c05c05e2ac"
    assert ledger["integrated_prs"] == [183, 184, 185, 186]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {183, 184, 185, 186}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "pressure chain" in rows[183]["capability"]
    assert "naturalRemainder(x0)" in rows[183]["remaining_boundary"]
    assert "k=1" in rows[184]["capability"]
    assert "Picard convergence" in rows[184]["remaining_boundary"]
    assert "exactly one" in rows[185]["capability"]
    assert "CellIndex" in rows[185]["remaining_boundary"]
    assert "multi-index Leibniz rule" in rows[186]["capability"]
    assert "provider inputs" in rows[186]["remaining_boundary"]


def test_source_addenda_and_runtime_stay_fail_closed() -> None:
    wide_pressure = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_wide_natural_pressure.json"
    )
    assert wide_pressure["full_reconstruction"] is False
    assert wide_pressure["paper_exact_velocity_available"] is False
    assert wide_pressure["layer"]["status"] == "formal-structure"
    assert wide_pressure["layer"]["actual_data"]["source_formula"] == "a^2 * phi0^2"
    assert "naturalRemainder(x0)" in wide_pressure["layer"]["remaining_boundary"]

    k1 = _load(ROOT / "references" / "provenance_manifest_addendum_184.json")
    assert k1["full_reconstruction"] is False
    assert k1["paper_exact_velocity_available"] is False
    assert k1["layers"][0]["status"] == "formal-structure"
    assert "W_n^(1)" in k1["layers"][0]["landed_capability"]
    assert "not Picard convergence" in k1["layers"][0]["remaining_boundary"]

    physical_scope = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_PHYSICAL_BASE_COVERAGE_PROVENANCE.json"
    )
    assert physical_scope["full_reconstruction"] is False
    assert physical_scope["paper_exact_velocity_available"] is False
    assert physical_scope["layer"]["status"] == "formal-structure"
    truth = physical_scope["layer"]["truth_boundary"]
    assert truth["declared_scope_physical_base_coverage_machine_checked"] is True
    assert truth["actual_cell_index_enumeration_machine_verified"] is False
    assert truth["actual_physical_base_values_materialized"] is False
    assert truth["actual_base_fields_verified"] is False
    assert truth["uniform_eq_7_9_to_7_11_verified"] is False
    assert truth["paper_exact_velocity_available"] is False

    weighted = _load(ROOT / "references" / "provenance_manifest_addendum_186.json")
    assert weighted["full_reconstruction"] is False
    assert weighted["paper_exact_velocity_available"] is False
    assert weighted["layers"][0]["status"] == "formal-structure"
    assert "Leibniz rule" in weighted["layers"][0]["capability"]
    assert "provider inputs" in weighted["layers"][0]["remaining_boundary"]
    assert any("infinite Eq. (9.21)" in item for item in weighted["forbidden_inferences"])

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_are_explicit() -> None:
    ledger = _load(LEDGER)
    forbidden = ledger["forbidden_inferences"]

    assert any("pressure branch" in item and "naturalRemainder(x0)" in item and "fixed point" in item for item in forbidden)
    assert any("k=1" in item and "Picard convergence" in item and "Proposition 5.3" in item for item in forbidden)
    assert any("CellIndex" in item and "FinalSlowBase.velocity" in item for item in forbidden)
    assert any("provider-supplied" in item and "infinite Eq. (9.21)" in item and "residual artifact" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)
