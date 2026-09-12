import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_178_181.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_178_181_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "527d2eb0989ea6104b5359c0ecd5330ff253d717"
    assert ledger["integrated_prs"] == [178, 179, 180, 181]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {178, 179, 180, 181}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "a^2*phi0^2" in rows[178]["capability"]
    assert "naturalRemainder(x0)" in rows[178]["remaining_boundary"]
    assert "k=0" in rows[179]["capability"]
    assert "k=1" in rows[179]["remaining_boundary"]
    assert "FinalSlowBase.velocity" in rows[180]["remaining_boundary"]
    assert "CellIndex" in rows[180]["remaining_boundary"]
    assert "provider-supplied" in rows[181]["capability"]
    assert "genuine Navier-Stokes residual artifact" in rows[181]["remaining_boundary"]


def test_source_addenda_and_runtime_remain_fail_closed() -> None:
    wide_source = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_wide_natural_source.json"
    )
    assert wide_source["full_reconstruction"] is False
    assert wide_source["paper_exact_velocity_available"] is False
    assert wide_source["layer"]["status"] == "formal-structure"
    assert wide_source["layer"]["actual_data"]["source_formula"] == "a^2 * phi0^2"
    assert "naturalRemainder(x0)" in wide_source["layer"]["remaining_boundary"]

    first_picard = _load(ROOT / "references" / "provenance_manifest_addendum_179.json")
    assert first_picard["full_reconstruction"] is False
    assert first_picard["paper_exact_velocity_available"] is False
    assert first_picard["layers"][0]["status"] == "formal-structure"
    assert "W_n^(0)=G f_n" in first_picard["layers"][0]["landed_capability"]
    assert "nontrivial k=1" in first_picard["layers"][0]["remaining_boundary"]

    physical_base = _load(
        ROOT / "references" / "provenance_manifest_addendum_phase_large_band_physical_base.json"
    )
    assert physical_base["full_reconstruction"] is False
    assert physical_base["paper_exact_velocity_available"] is False
    assert physical_base["layer"]["status"] == "formal-structure"
    assert "FinalSlowBase.velocity" in physical_base["layer"]["capability"]
    assert "does not evaluate the noncomputable Prepared object" in physical_base["layer"]["remaining_boundary"]

    prefix_jet = _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_eq921_prefix_jet.json"
    )
    assert prefix_jet["full_reconstruction"] is False
    assert prefix_jet["paper_exact_velocity_available"] is False
    assert prefix_jet["layer"]["status"] == "formal-structure"
    flags = prefix_jet["truth_flags"]
    assert flags["analytic_derivatives_machine_derived_from_actual_corrections"] is False
    assert flags["paper_fixed_cutoff_derivatives_machine_verified"] is False
    assert flags["actual_section9_sequence_verified"] is False
    assert flags["eq_9_21_infinite_sum_constructed"] is False
    assert flags["endpoint_covered"] is False
    assert flags["residual_artifact_ready"] is False
    assert flags["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_bindings_preserve_independent_verification_boundaries() -> None:
    ledger = _load(LEDGER)
    rows = {row["pr"]: row for row in ledger["layers"]}

    assert "not a complete naturalRemainder(x0)" in rows[178]["remaining_boundary"]
    assert "Picard convergence" in rows[179]["remaining_boundary"]
    assert "does not evaluate the noncomputable Prepared object" in rows[180]["remaining_boundary"]
    assert "not machine-derived" in rows[181]["remaining_boundary"]
    assert "infinite Eq. (9.21)" in rows[181]["remaining_boundary"]
    assert "smooth compact forcing" in rows[181]["remaining_boundary"]

    forbidden = ledger["forbidden_inferences"]
    assert any("naturalRemainder(x0)" in item and "fixed point" in item for item in forbidden)
    assert any("k=0" in item and "k=1" in item and "Proposition 5.3" in item for item in forbidden)
    assert any("Prepared" in item and "CellIndex" in item for item in forbidden)
    assert any("provider-supplied" in item and "residual artifact" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)
