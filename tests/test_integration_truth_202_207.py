import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_202_207.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_202_207_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "e4b17581a915e230e506f6beeb5cd19789fdea37"
    assert ledger["integrated_prs"] == [202, 203, 204, 205, 206, 207]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {202, 203, 204, 205, 206, 207}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "first genuine theorem-scale Picard iterate x1" in rows[202]["capability"]
    assert "naturalRemainder(x1)" in rows[202]["remaining_boundary"]
    assert "content-addressed fail-closed handoff" in rows[203]["capability"]
    assert "CellIndex" in rows[203]["remaining_boundary"]
    assert "timeSwitch" in rows[204]["equation_or_source"]
    assert "smooth field extension through t=1" in rows[204]["remaining_boundary"]
    assert "1/(n+1)" in rows[205]["capability"]
    assert "naturalRemainder(x1)" in rows[205]["remaining_boundary"]
    assert "preceding-diffusion" in rows[206]["capability"]
    assert "hierarchy-owned partial_eta^2 actualLowerSource" in rows[206]["remaining_boundary"]
    assert "common_curl_and_divergence" in rows[207]["equation_or_source"]
    assert "not machine replayed" in rows[207]["remaining_boundary"]


def test_landed_pr_202_207_truth_boundaries_remain_fail_closed() -> None:
    first_picard = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard.json"
    )
    assert first_picard["full_reconstruction"] is False
    assert first_picard["paper_exact_velocity_available"] is False
    assert first_picard["layer"]["status"] == "formal-structure"
    assert "first genuine theorem-scale Picard x1" in first_picard["layer"]["capability"]
    assert "fixed-point" in first_picard["layer"]["remaining_boundary"]

    canonical_export = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_CANONICAL_APPLICATION_EXPORT_PROVENANCE.json"
    )
    assert canonical_export["full_reconstruction"] is False
    assert canonical_export["paper_exact_velocity_available"] is False
    truth203 = canonical_export["layer"]["truth_boundary"]
    assert truth203["actual_canonical_prepared_application_machine_verified"] is False
    assert truth203["actual_cell_index_enumeration_machine_verified"] is False
    assert truth203["actual_physical_base_values_materialized"] is False
    assert truth203["theorem_applications_machine_replayed"] is False
    assert truth203["uniform_eq_7_9_to_7_11_verified"] is False

    endpoint_switch = (
        ROOT / "references" / "SECTION10_PAPER_TIME_SWITCH_ENDPOINT_JET_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "no actual Section 9 all-order one-sided endpoint jet" in endpoint_switch
    assert "no smooth extension through `t=1`" in endpoint_switch
    assert "endpoint_residual_closure_verified" in endpoint_switch
    assert "paper_exact_velocity_available" in endpoint_switch

    first_picard_average = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_average.json"
    )
    assert first_picard_average["full_reconstruction"] is False
    assert first_picard_average["paper_exact_velocity_available"] is False
    assert first_picard_average["layer"]["status"] == "formal-structure"
    assert "coefficientOperators.average(u1)" in first_picard_average["layer"]["capability"]
    assert "naturalRemainder(x1)" in first_picard_average["layer"]["remaining_boundary"]

    preceding_diffusion = _load(ROOT / "references" / "provenance_manifest_addendum_206.json")
    assert preceding_diffusion["full_reconstruction"] is False
    assert preceding_diffusion["paper_exact_velocity_available"] is False
    assert preceding_diffusion["layers"][0]["status"] == "formal-structure"
    assert "fourth-mixed" in preceding_diffusion["scope_note"]
    assert "hierarchy-owned partial_eta^2 actualLowerSource" in preceding_diffusion["layers"][0]["remaining_boundary"]

    common_curl = _load(ROOT / "references" / "ACTUAL_SIGNED_COMMON_CURL_PROVENANCE.json")
    assert common_curl["status"] == "formal-structure"
    truth207 = common_curl["truth_boundary"]
    assert truth207["theorem_application_machine_replayed"] is False
    assert truth207["actual_signed_wave_values_materialized"] is False
    assert truth207["canonical_scope_link_verified"] is False
    assert truth207["amplitude_ode_inputs_verified"] is False
    assert truth207["paper_exact_divergence_free_wave_available"] is False
    assert truth207["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_are_explicit() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any(
        "first Picard x1" in item
        and "naturalRemainder(x1)" in item
        and "Picard convergence" in item
        for item in forbidden
    )
    assert any(
        "coefficientOperators.average(u1)" in item
        and "NaturalProfileAssembly" in item
        for item in forbidden
    )
    assert any(
        "canonical Prepared export" in item
        and "CellIndex" in item
        and "Eqs. (7.9)-(7.11)" in item
        for item in forbidden
    )
    assert any(
        "timeSwitch" in item
        and "Section 9" in item
        and "final forcing" in item
        for item in forbidden
    )
    assert any(
        "preceding-diffusion" in item
        and "hierarchy ownership" in item
        and "Proposition 5.3" in item
        for item in forbidden
    )
    assert any(
        "common_curl_and_divergence" in item
        and "Lean theorem" in item
        and "paper-exact" in item
        for item in forbidden
    )
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
