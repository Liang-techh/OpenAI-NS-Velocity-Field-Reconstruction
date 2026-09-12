import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_220_221.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_220_221_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "4b01ec10be9eed10d6c1737f5f2a7f237c9e4912"
    assert ledger["integrated_prs"] == [220, 221]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {220, 221}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()
        if row["provenance"] is not None:
            assert (ROOT / row["provenance"]).is_file()

    assert "Prepared.N + 1 <= n" in rows[220]["capability"]
    assert "finite-head cutoff deficit" in rows[220]["remaining_boundary"]
    assert "Cartesian space" in rows[221]["capability"]
    assert "fails closed on the symmetry axis" in rows[221]["capability"]
    assert rows[221]["provenance"] is None
    assert "truth boundary" in rows[221]["provenance_note"]
    assert "provider-supplied" in rows[221]["remaining_boundary"]
    assert "genuine Navier-Stokes residual/forcing" in rows[221]["remaining_boundary"]


def test_landed_truth_sources_for_220_221_remain_fail_closed() -> None:
    mean_tail = _load(ROOT / "references" / "ACTUAL_SIGNED_MEAN_TAIL_PROVENANCE.json")
    assert mean_tail["full_reconstruction"] is False
    assert mean_tail["paper_exact_velocity_available"] is False
    assert mean_tail["layer"]["status"] == "formal-structure"
    truth220 = mean_tail["layer"]["truth_boundary"]
    assert truth220["literal_requested_cross_tail_theorem_machine_replayed"] is False
    assert truth220["actual_mean_cross_values_materialized"] is False
    assert truth220["finite_head_mean_defect_solved"] is False
    assert truth220["compact_mean_correction_available"] is False
    assert truth220["paper_exact_velocity_available"] is False

    jacobian_source = (
        ROOT
        / "src"
        / "openai_ns_reconstruction"
        / "section9_section10_localized_velocity_spatial_jacobian.py"
    ).read_text(encoding="utf-8")
    assert "axis_spatial_derivative_covered: bool = False" in jacobian_source
    assert "actual_section9_sequence_verified: bool = False" in jacobian_source
    assert "eq_9_21_infinite_sum_constructed: bool = False" in jacobian_source
    assert "section9_field_smooth_extension_through_t1_constructed: bool = False" in jacobian_source
    assert "residual_artifact_ready: bool = False" in jacobian_source
    assert "endpoint_residual_closure_verified: bool = False" in jacobian_source
    assert "paper_exact_velocity_available: bool = False" in jacobian_source
    assert "Production code performs no finite differencing" in jacobian_source
    assert "axis regularity must be supplied separately" in jacobian_source

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
    assert any("mean-cross identity" in x and "finite-head mean defect" in x for x in forbidden)
    assert any("off-axis spatial Jacobian" in x and "Navier-Stokes residual or forcing" in x for x in forbidden)
    assert any("symmetry axis" in x and "axis-regularity" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
