import json
from pathlib import Path

from openai_ns_reconstruction.provenance import (
    IncompleteReconstructionError,
    require_complete_reconstruction,
)
from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_235_236.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_235_236_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "c7fabfe1798ec8e416c543556bf91348a7bceebe"
    assert ledger["integrated_prs"] == [235, 236]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {235, 236}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "ProfileFourthMixedJet" in rows[235]["capability"]
    assert "partial_eta^2 actualLowerSource" in rows[235]["remaining_boundary"]
    assert "rank_rows_on_patch" in rows[236]["capability"]
    assert "CorrectionState.debt" in rows[236]["capability"]
    assert "finite-head signed-defect bridge" in rows[236]["remaining_boundary"]


def test_pr236_source_truth_remains_fail_closed() -> None:
    source = _load(ROOT / "references" / "COMPACT_MEAN_CORRECTION_FORMAL_PROVENANCE.json")
    assert source["full_reconstruction"] is False
    assert source["paper_exact_velocity_available"] is False
    layer = source["layer"]
    assert layer["status"] == "formal-structure"
    truth = layer["truth_boundary"]
    assert truth["theorem_application_machine_replayed"] is False
    assert truth["actual_debt_values_materialized"] is False
    assert truth["actual_bump_values_materialized"] is False
    assert truth["finite_head_signed_defect_link_verified"] is False
    assert truth["compact_mean_correction_field_materialized"] is False
    assert truth["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"

    try:
        require_complete_reconstruction()
    except IncompleteReconstructionError:
        pass
    else:
        raise AssertionError("paper-exact completion gate must remain fail-closed")


def test_forbidden_inferences_and_agent_contract_stay_locked() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("partial_eta^2 actualLowerSource/forcing" in item for item in forbidden)
    assert any("rank_rows_on_patch" in item and "Mathlib bump" in item for item in forbidden)
    assert any("signed mean defect" in item and "CorrectionState.debt" in item for item in forbidden)
    assert any("surrogate bump" in item and "paper-exact" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
