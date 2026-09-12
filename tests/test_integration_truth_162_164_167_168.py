import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_162_164_167_168.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_latest_stage1_3_increments_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "3dbe9002c440ee01cc6eb72523c842b905dd8fb3"
    assert ledger["integrated_prs"] == [162, 163, 164, 167, 168]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {162, 163, 164, 167, 168}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "naturalRemainder(x0)" in rows[162]["remaining_boundary"]
    assert "genuine Picard x1/fixed point" in rows[162]["remaining_boundary"]
    assert "analytic eta jet of final Eq. (5.7) forcing f_n" in rows[163]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in rows[163]["remaining_boundary"]
    assert "Proposition 5.5 paper-exact base provider" in rows[164]["remaining_boundary"]
    assert "curl-realized oscillatory wave" in rows[164]["remaining_boundary"]
    assert "theorem-selected t=1/Lambda" in rows[167]["remaining_boundary"]
    assert "naturalRemainder(x0)" in rows[167]["remaining_boundary"]
    assert "analytic eta jet of the final Eq. (5.7) forcing f_n" in rows[168]["remaining_boundary"]
    assert "genuine hierarchy-owned k=1 Picard application" in rows[168]["remaining_boundary"]


def test_source_addenda_canonical_manifest_and_runtime_stay_fail_closed() -> None:
    stage1_resolvent = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_natural_resolvent.json"
    )
    assert stage1_resolvent["full_reconstruction"] is False
    assert stage1_resolvent["paper_exact_velocity_available"] is False
    assert stage1_resolvent["layer"]["status"] == "formal-structure"

    stage2_history = _load(ROOT / "references" / "provenance_manifest_addendum_163.json")
    assert stage2_history["full_reconstruction"] is False
    assert stage2_history["paper_exact_velocity_available"] is False
    assert stage2_history["layers"][0]["status"] == "formal-structure"

    stage3 = _load(
        ROOT / "references" / "provenance_manifest_addendum_phase_large_band_family_coverage.json"
    )
    assert stage3["full_reconstruction"] is False
    assert stage3["paper_exact_velocity_available"] is False
    assert stage3["layer"]["status"] == "formal-structure"

    stage1_remainder = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_natural_remainder.json"
    )
    assert stage1_remainder["full_reconstruction"] is False
    assert stage1_remainder["paper_exact_velocity_available"] is False
    assert stage1_remainder["layer"]["status"] == "formal-structure"
    actual = stage1_remainder["layer"]["actual_data"]
    assert actual["coefficient_operators_caller_supplied"] is False
    assert actual["axis_data_caller_supplied"] is False
    assert actual["natural_resolvent_caller_supplied"] is False
    assert actual["t_and_a_are_theorem_arguments_not_selected_here"] is True

    stage2_parameter = _load(ROOT / "references" / "provenance_manifest_addendum_168.json")
    assert stage2_parameter["full_reconstruction"] is False
    assert stage2_parameter["paper_exact_velocity_available"] is False
    assert stage2_parameter["layers"][0]["status"] == "formal-structure"
    assert "final Eq. (5.7) forcing f_n" in stage2_parameter["layers"][0]["remaining_boundary"]

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[3]["status"] == "formal-structure"


def test_new_composition_and_parameter_jet_do_not_close_picard_or_wave_gates() -> None:
    ledger = _load(LEDGER)
    rows = {row["pr"]: row for row in ledger["layers"]}

    assert "exact finite alternating prefix" in rows[162]["capability"]
    assert "pinned naturalRemainder expression is now executable" in rows[167]["capability"]
    assert "normalized amplitude coefficient state" in rows[167]["remaining_boundary"]
    assert "genuine first Picard iterate/fixed point" in rows[167]["remaining_boundary"]
    assert "global all-index weighted AxisSpace membership/norms" in rows[167]["remaining_boundary"]

    assert "analytic eta derivative" in rows[168]["capability"]
    assert "final Eq. (5.7) forcing f_n" in rows[168]["remaining_boundary"]
    assert "recursive coefficient materialization" in rows[168]["remaining_boundary"]

    forbidden = ledger["forbidden_inferences"]
    assert any("naturalRemainder(x0)" in item and "fixed point" in item for item in forbidden)
    assert any("strict-lower-source eta jet" in item and "final f_n" in item for item in forbidden)
    assert any("finite active-family admission" in item for item in forbidden)
    assert any("setting f equal to the evaluated residual" in item for item in forbidden)
