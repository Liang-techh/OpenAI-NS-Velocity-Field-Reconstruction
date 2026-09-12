import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_173_176.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_173_176_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "d9a9bf88f8f8a75c43c7a32130a3de776f0fbf98"
    assert ledger["integrated_prs"] == [173, 174, 175, 176]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {173, 174, 175, 176}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "naturalRemainder(x0)" in rows[173]["remaining_boundary"]
    assert "A0 W + A1 partial_eta W" in rows[174]["remaining_boundary"]
    assert "CellIndex" in rows[175]["remaining_boundary"]
    assert "source-not-residual" in rows[176]["capability"]
    assert "not the genuine Eq. (9.21) residual artifact" in rows[176]["remaining_boundary"]


def test_source_addenda_and_runtime_remain_fail_closed() -> None:
    amplitude = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_amplitude.json"
    )
    assert amplitude["full_reconstruction"] is False
    assert amplitude["paper_exact_velocity_available"] is False
    assert amplitude["layer"]["status"] == "formal-structure"
    assert amplitude["layer"]["actual_data"]["magnitude_representation"] == "signed log"
    assert "naturalRemainder(x0)" in amplitude["layer"]["remaining_boundary"]

    forcing = _load(ROOT / "references" / "provenance_manifest_addendum_174.json")
    assert forcing["full_reconstruction"] is False
    assert forcing["paper_exact_velocity_available"] is False
    assert forcing["layers"][0]["status"] == "formal-structure"
    assert "partial_eta f_n" in forcing["layers"][0]["landed_capability"]

    prepared = _load(
        ROOT / "references" / "provenance_manifest_addendum_phase_large_band_prepared_source.json"
    )
    assert prepared["full_reconstruction"] is False
    assert prepared["paper_exact_velocity_available"] is False
    assert prepared["layer"]["status"] == "formal-structure"
    assert "does not machine-discharge exists_prepared hypotheses" in prepared["layer"]["remaining_boundary"]

    prefix = _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_eq921_prefix_source_artifact.json"
    )
    assert prefix["status"] == "formal-structure"
    assert prefix["finite_prefix_source_materialized"] is True
    assert prefix["residual_artifact_ready"] is False
    assert prefix["source_majorants_derived_from_actual_residual_verified"] is False
    assert prefix["actual_section9_sequence_verified"] is False
    assert prefix["eq_9_21_infinite_sum_constructed"] is False
    assert prefix["endpoint_limits_constructed"] is False
    assert prefix["all_order_borel_smoothness_verified"] is False
    assert prefix["smooth_compact_forcing_verified"] is False
    assert prefix["full_reconstruction"] is False
    assert prefix["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_latest_bindings_do_not_close_fixed_point_hierarchy_phase_or_endpoint_gates() -> None:
    ledger = _load(LEDGER)
    rows = {row["pr"]: row for row in ledger["layers"]}

    assert "signed-log" in rows[173]["capability"]
    assert "fixed point" in rows[173]["remaining_boundary"]
    assert "partial_eta f_n" in rows[174]["capability"]
    assert "Proposition 5.3" in rows[174]["remaining_boundary"]
    assert "Prepared object" in rows[175]["remaining_boundary"]
    assert "Proposition 5.5" in rows[175]["remaining_boundary"]
    assert "externally supplied non-paper correction values" in rows[176]["remaining_boundary"]
    assert "all-order Borel smoothness" in rows[176]["remaining_boundary"]
    assert "smooth compact forcing" in rows[176]["remaining_boundary"]

    forbidden = ledger["forbidden_inferences"]
    assert any("naturalRemainder(x0)" in item for item in forbidden)
    assert any("partial_eta f_n" in item and "Proposition 5.3" in item for item in forbidden)
    assert any("Prepared object" in item and "CellIndex" in item for item in forbidden)
    assert any("source-prefix bytes" in item and "residual artifact" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)
