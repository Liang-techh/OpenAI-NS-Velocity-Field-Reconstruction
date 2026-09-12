import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_215_218.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_215_218_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "0b2d0b56b350939196cfe83defa0a9afeaca82c4"
    assert ledger["integrated_prs"] == [215, 216, 217, 218]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {215, 216, 217, 218}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()
        if row["provenance"] is not None:
            assert (ROOT / row["provenance"]).is_file()

    assert "u1 * u1" in rows[215]["capability"]
    assert "naturalRemainder(x1)" in rows[215]["remaining_boundary"]
    assert "second eta derivative" in rows[216]["capability"]
    assert "fourth-mixed" in rows[216]["remaining_boundary"]
    assert "canonical Prepared signed-label roster" in rows[217]["capability"]
    assert "not machine replayed" in rows[217]["remaining_boundary"]
    assert "analytic Cartesian Hessian" in rows[218]["capability"]
    assert rows[218]["provenance"] is None
    assert "machine-readable truth record" in rows[218]["provenance_note"]
    assert "transition collar" in rows[218]["remaining_boundary"]


def test_landed_truth_sources_remain_fail_closed() -> None:
    axial_square = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_axial_square.json"
    )
    assert axial_square["full_reconstruction"] is False
    assert axial_square["paper_exact_velocity_available"] is False
    assert axial_square["layer"]["status"] == "formal-structure"
    assert "naturalRemainder(x1)" in axial_square["layer"]["remaining_boundary"]

    omega_second = _load(ROOT / "references" / "provenance_manifest_addendum_216.json")
    assert omega_second["full_reconstruction"] is False
    assert omega_second["paper_exact_velocity_available"] is False
    assert omega_second["layers"][0]["status"] == "formal-structure"
    assert "hierarchy-owned partial_eta^2 actualLowerSource" in omega_second["layers"][0]["remaining_boundary"]

    signed_scope = _load(ROOT / "references" / "ACTUAL_SIGNED_CANONICAL_SCOPE_PROVENANCE.json")
    truth217 = signed_scope["layer"]["truth_boundary"]
    assert truth217["actual_signed_label_export_machine_verified"] is False
    assert truth217["canonical_prepared_application_machine_replayed"] is False
    assert truth217["actual_cycle_invariant_machine_verified"] is False
    assert truth217["common_curl_theorem_machine_replayed"] is False
    assert truth217["actual_signed_wave_values_materialized"] is False
    assert truth217["paper_exact_divergence_free_wave_available"] is False
    assert truth217["paper_exact_velocity_available"] is False

    localization_source = (
        ROOT / "src" / "openai_ns_reconstruction" / "spatial_localization.py"
    ).read_text(encoding="utf-8")
    assert "section10_spatial_cutoff_hessian" in localization_source
    assert "not promote the executable transition collar" in localization_source

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_and_agent_contract() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("u1 * u1" in x and "naturalRemainder(x1)" in x for x in forbidden)
    assert any("second-eta" in x and "Proposition 5.3" in x for x in forbidden)
    assert any("canonical Prepared signed-label roster" in x and "replay Lean" in x for x in forbidden)
    assert any("analytic Hessian" in x and "Navier-Stokes residual" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
