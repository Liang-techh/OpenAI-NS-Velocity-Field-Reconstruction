import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_248_255.json"
PR255_PROVENANCE = ROOT / "references" / "provenance_manifest_addendum_255.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_records_current_landed_truth_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "05ef960f29d8b9b21955c0eacbe9132a610c4cb1"
    assert ledger["integrated_prs"] == [248, 249, 250, 251, 253, 255]
    assert ledger["excluded_nonlanded_prs"] == [252, 254, 256]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {248, 249, 250, 251, 253, 255}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["landed_capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "dot1" in rows[248]["landed_capability"]
    assert "partial_eta^2 f_n" in rows[249]["landed_capability"]
    assert "LocalResidualFlatness.selected_schedule" in rows[250]["landed_capability"]
    assert "convective term" in rows[251]["landed_capability"]
    assert "mixed1(average(u1),u1)" in rows[253]["landed_capability"]
    assert "formal_operator_sum" in rows[255]["landed_capability"]
    assert "not a residual or forcing artifact" in rows[255]["remaining_boundary"]


def test_pr255_formal_operator_provenance_matches_fail_closed_source_contract() -> None:
    provenance = _load(PR255_PROVENANCE)
    assert provenance["source_pr"] == 255
    assert provenance["full_reconstruction"] is False
    assert provenance["paper_exact_velocity_available"] is False

    layer = provenance["layer"]
    assert layer["status"] == "formal-structure"
    boundary = layer["truth_boundary"]
    assert boundary["provider_supplied_finite_eq921_prefix"] is True
    assert boundary["off_axis_only"] is True
    assert boundary["formal_ns_operator_combination_verified"] is True
    for key in (
        "actual_section7_8_correction_sequence_constructed",
        "eq_9_21_infinite_sum_constructed",
        "axis_regularity_from_actual_paper_field_constructed",
        "section9_field_smooth_extension_through_t1_constructed",
        "residual_artifact_ready",
        "forcing_artifact_ready",
        "endpoint_residual_closure_verified",
        "finite_energy_closure_verified",
        "blow_up_closure_verified",
        "paper_exact_velocity_available",
    ):
        assert boundary[key] is False

    source = (ROOT / layer["artifact"]).read_text(encoding="utf-8")
    assert "formal_operator_sum" in source
    assert "residual_artifact_ready: bool = False" in source
    assert "forcing_artifact_ready: bool = False" in source
    assert "paper_exact_velocity_available: bool = False" in source
    assert "does not establish" in source or "does not" in source


def test_canonical_runtime_and_agent_contract_stay_fail_closed() -> None:
    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[6]["status"] == "formal-structure"

    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("slow2(x1)" in item and "fixed-point convergence" in item for item in forbidden)
    assert any("partial_eta^2 f_n" in item and "Proposition 5.3" in item for item in forbidden)
    assert any("selected_schedule" in item and "actual stage fields" in item for item in forbidden)
    assert any("formal_operator_sum" in item and "not a genuine residual" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
