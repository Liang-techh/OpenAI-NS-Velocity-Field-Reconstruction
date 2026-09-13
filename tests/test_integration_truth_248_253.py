import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_248_253.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_records_only_landed_prs_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "23025263693addd44f6443dfe95e81bd5940c535"
    assert ledger["integrated_prs"] == [248, 249, 250, 251, 253]
    assert ledger["excluded_nonlanded_prs"] == [252]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {248, 249, 250, 251, 253}
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
    assert "-param1(u1, AxisData.d*u1)" in rows[253]["remaining_boundary"]
    assert "slow2(x1)" in rows[253]["remaining_boundary"]
    assert "naturalRemainder(x1)" in rows[253]["remaining_boundary"]


def test_stage1_new_branch_and_canonical_runtime_stay_fail_closed() -> None:
    stage1 = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_average_mixed.json"
    )
    assert stage1["full_reconstruction"] is False
    assert stage1["paper_exact_velocity_available"] is False
    layer = stage1["layer"]
    assert layer["status"] == "formal-structure"
    boundary = layer["truth_boundary"]
    assert boundary["first_picard_x1_materialized"] is True
    assert boundary["slow2_average_mixed_branch_materialized"] is True
    assert boundary["slow2_param_branch_materialized"] is False
    assert boundary["slow2_complete"] is False
    assert boundary["natural_remainder_x1_materialized"] is False
    assert boundary["picard_x2_materialized"] is False
    assert boundary["fixed_point_materialized"] is False
    assert boundary["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[6]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_remain_explicit() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("slow2(x1)" in x and "fixed-point convergence" in x for x in forbidden)
    assert any("partial_eta^2 f_n" in x and "Proposition 5.3" in x for x in forbidden)
    assert any("selected_schedule" in x and "actual stage fields" in x for x in forbidden)
    assert any("convective term" in x and "full Navier-Stokes residual" in x for x in forbidden)
    assert any("PR #252" in x and "must not be counted" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
