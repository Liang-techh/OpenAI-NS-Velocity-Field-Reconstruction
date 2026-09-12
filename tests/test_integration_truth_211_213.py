import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_211_213.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_211_213_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "e96b6d965115929a3f99cda70d1473803f15bf34"
    assert ledger["integrated_prs"] == [211, 212, 213]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {211, 212, 213}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "phi1 * phi1" in rows[211]["capability"]
    assert "naturalRemainder(x1)" in rows[211]["remaining_boundary"]
    assert "fullRequest_jets_from_invariant" in rows[212]["capability"]
    assert "No Lean theorem is replayed" in rows[212]["remaining_boundary"]
    assert "time derivative d_t u_cut" in rows[213]["capability"]
    assert "Navier-Stokes residual artifact" in rows[213]["remaining_boundary"]


def test_landed_truth_sources_remain_fail_closed() -> None:
    square = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_angular_square.json"
    )
    assert square["full_reconstruction"] is False
    assert square["paper_exact_velocity_available"] is False
    assert square["layer"]["status"] == "formal-structure"
    assert "naturalRemainder(x1)" in square["layer"]["remaining_boundary"]

    invariant = _load(
        ROOT / "references" / "provenance_manifest_addendum_actual_signed_common_curl_invariant.json"
    )
    truth212 = invariant["layer"]["truth_boundary"]
    assert truth212["full_request_jets_theorem_machine_replayed"] is False
    assert truth212["actual_cycle_invariant_machine_verified"] is False
    assert truth212["common_curl_theorem_machine_replayed"] is False
    assert truth212["actual_signed_wave_values_materialized"] is False
    assert truth212["paper_exact_divergence_free_wave_available"] is False
    assert truth212["paper_exact_velocity_available"] is False

    time_derivative = _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_section10_time_derivative.json"
    )
    truth213 = time_derivative["layer"]["truth_boundary"]
    assert truth213["provider_supplied_finite_prefix_jets"] is True
    assert truth213["production_finite_difference_used"] is False
    assert truth213["actual_section9_sequence_verified"] is False
    assert truth213["eq_9_21_infinite_sum_constructed"] is False
    assert truth213["section9_field_smooth_extension_through_t1_constructed"] is False
    assert truth213["residual_artifact_ready"] is False
    assert truth213["endpoint_residual_closure_verified"] is False
    assert truth213["paper_exact_velocity_available"] is False

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
    assert any("phi1 * phi1" in x and "Picard convergence" in x for x in forbidden)
    assert any("fullRequest_jets_from_invariant" in x and "replay Lean" in x for x in forbidden)
    assert any("finite provider-supplied localized prefix" in x and "residual/forcing" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
