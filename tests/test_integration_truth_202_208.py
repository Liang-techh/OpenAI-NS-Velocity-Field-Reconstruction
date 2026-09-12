import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_202_208.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_202_208_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "2d74857236cac1f3e553e82c040e1bee386e2d89"
    assert ledger["integrated_prs"] == [202, 203, 204, 205, 206, 207, 208]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {202, 203, 204, 205, 206, 207, 208}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "first genuine theorem-scale Picard iterate x1" in rows[202]["capability"]
    assert "naturalRemainder(x1)" in rows[202]["remaining_boundary"]
    assert "CellIndex" in rows[203]["remaining_boundary"]
    assert "timeSwitch" in rows[204]["capability"]
    assert "smooth extension through t=1" in rows[204]["remaining_boundary"]
    assert "1/(n+1)" in rows[205]["capability"]
    assert "partial_eta^2 actualLowerSource" in rows[206]["remaining_boundary"]
    assert "common_curl_and_divergence" in rows[207]["capability"]
    assert "not machine replayed" in rows[207]["remaining_boundary"]
    assert "u_cut = c curl(A) + grad(c) x A + c B e_theta" in rows[208]["capability"]
    assert "residual/forcing closure" in rows[208]["remaining_boundary"]


def test_landed_truth_sources_remain_fail_closed() -> None:
    first_picard = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_wide_first_picard.json"
    )
    assert first_picard["full_reconstruction"] is False
    assert first_picard["paper_exact_velocity_available"] is False
    assert first_picard["layer"]["status"] == "formal-structure"

    first_average = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_average.json"
    )
    assert first_average["paper_exact_velocity_available"] is False
    assert "naturalRemainder(x1)" in first_average["layer"]["remaining_boundary"]

    canonical_export = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_CANONICAL_APPLICATION_EXPORT_PROVENANCE.json"
    )
    truth203 = canonical_export["layer"]["truth_boundary"]
    assert truth203["actual_canonical_prepared_application_machine_verified"] is False
    assert truth203["actual_cell_index_enumeration_machine_verified"] is False
    assert truth203["actual_physical_base_values_materialized"] is False
    assert truth203["theorem_applications_machine_replayed"] is False
    assert truth203["paper_exact_velocity_available"] is False

    preceding = _load(ROOT / "references" / "provenance_manifest_addendum_206.json")
    assert preceding["full_reconstruction"] is False
    assert preceding["paper_exact_velocity_available"] is False
    assert "hierarchy-owned partial_eta^2 actualLowerSource" in preceding["layers"][0]["remaining_boundary"]

    common_curl = _load(ROOT / "references" / "ACTUAL_SIGNED_COMMON_CURL_PROVENANCE.json")
    truth207 = common_curl["truth_boundary"]
    assert truth207["theorem_application_machine_replayed"] is False
    assert truth207["actual_signed_wave_values_materialized"] is False
    assert truth207["paper_exact_divergence_free_wave_available"] is False
    assert truth207["paper_exact_velocity_available"] is False

    endpoint_switch = (
        ROOT / "references" / "SECTION10_PAPER_TIME_SWITCH_ENDPOINT_JET_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "no actual Section 9 all-order one-sided endpoint jet" in endpoint_switch
    assert "no smooth extension through `t=1`" in endpoint_switch

    localized = (
        ROOT / "references" / "SECTION9_SECTION10_LOCALIZED_VELOCITY_JET_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "finite-prefix analytic adapter only" in localized
    assert "does not construct a Navier--Stokes residual or forcing" in localized
    assert "point_binding_to_actual_correction_field_verified" in localized
    assert "actual_section9_sequence_verified" in localized
    assert "eq_9_21_infinite_sum_constructed" in localized
    assert "residual_artifact_ready" in localized
    assert "paper_exact_velocity_available" in localized

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
    assert any("first Picard x1" in x and "Picard convergence" in x for x in forbidden)
    assert any("canonical Prepared export" in x and "CellIndex" in x for x in forbidden)
    assert any("timeSwitch" in x and "final forcing" in x for x in forbidden)
    assert any("preceding-diffusion" in x and "Proposition 5.3" in x for x in forbidden)
    assert any("common_curl_and_divergence" in x and "materialize" in x for x in forbidden)
    assert any("finite Eq. (9.21) prefix jet" in x and "bounded energy" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
