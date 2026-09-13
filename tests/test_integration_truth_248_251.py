import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_248_251.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_248_251_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "ac9b2f065e2ad6be57a324941a840b8372348a96"
    assert ledger["integrated_prs"] == [248, 249, 250, 251]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False
    assert "squash commit subject for PR #250 is misleading" in ledger["historical_note"]
    assert "selected-schedule residual-flatness" in ledger["historical_note"]

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {248, 249, 250, 251}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["landed_capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "dot1((2*D*eta)*average(u1),u1)" in rows[248]["landed_capability"]
    assert "not complete slow2(x1) or naturalRemainder(x1)" in rows[248]["remaining_boundary"]
    assert "partial_eta^2 f_n" in rows[249]["landed_capability"]
    assert "Proposition 5.3 all-jets residual decay" in rows[249]["remaining_boundary"]
    assert "LocalResidualFlatness.selected_schedule" in rows[250]["landed_capability"]
    assert "not machine-materialized" in rows[250]["remaining_boundary"]
    assert "(u_cut dot grad)u_cut" in rows[251]["landed_capability"]
    assert "genuine Navier-Stokes residual" in rows[251]["remaining_boundary"]
    assert "paper-exact velocity" in rows[251]["remaining_boundary"]


def test_landed_provenance_for_248_251_remains_fail_closed() -> None:
    stage1 = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_average_dot.json"
    )
    assert stage1["full_reconstruction"] is False
    assert stage1["paper_exact_velocity_available"] is False
    assert stage1["layer"]["status"] == "formal-structure"
    assert "complete slow2(x1)" in stage1["layer"]["validation"][-1]
    assert "fixed-point convergence/finality" in stage1["layer"]["validation"][-1]

    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_249.json")
    assert stage2["full_reconstruction"] is False
    assert stage2["paper_exact_velocity_available"] is False
    assert stage2["layers"][0]["status"] == "formal-structure"
    assert "partial_eta W_n^(1)" in stage2["layers"][0]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in stage2["layers"][0]["remaining_boundary"]

    section9 = _load(
        ROOT / "references" / "SECTION9_SELECTED_SCHEDULE_FLATNESS_PROVENANCE.json"
    )
    assert section9["full_reconstruction"] is False
    assert section9["paper_exact_velocity_available"] is False
    assert section9["layer"]["status"] == "formal-structure"
    boundary = section9["layer"]["truth_boundary"]
    assert boundary["selected_schedule_theorem_machine_replayed"] is False
    assert boundary["selected_schedule_values_materialized"] is False
    assert boundary["finite_stage_fields_materialized"] is False
    assert boundary["residual_jet_rate_values_materialized"] is False
    assert boundary["section9_iteration_machine_materialized"] is False
    assert boundary["paper_exact_velocity_available"] is False

    convection = (
        ROOT / "references" / "SECTION10_CONVECTIVE_ACCELERATION_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "formal-structure analytic operator adapter only" in convection
    assert "provider-supplied" in convection
    assert "residual_artifact_ready=false" in convection
    assert "paper_exact_velocity_available=false" in convection
    assert "symmetry-axis regularity" in convection

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[6]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_and_agent_contract() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("slow2(x1)" in x and "fixed-point convergence" in x for x in forbidden)
    assert any("partial_eta^2 f_n" in x and "Proposition 5.3" in x for x in forbidden)
    assert any("selected_schedule" in x and "actual stage fields" in x for x in forbidden)
    assert any("convective term" in x and "full Navier-Stokes residual" in x for x in forbidden)
    assert any("PR #250" in x and "not provenance evidence" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
