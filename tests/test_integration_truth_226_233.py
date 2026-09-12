import json
from pathlib import Path

from openai_ns_reconstruction.provenance import (
    IncompleteReconstructionError,
    require_complete_reconstruction,
)
from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_226_233.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_landed_prs_226_233_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "a05356dd69869b997d45e0db8a3896ace3aa0faf"
    assert ledger["integrated_prs"] == [226, 227, 228, 229, 231, 232, 233]
    assert 230 not in ledger["integrated_prs"]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {226, 227, 228, 229, 231, 232, 233}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()
        if row["provenance"] is not None:
            assert (ROOT / row["provenance"]).is_file()

    assert "first pinned slow2(x1) branch" in rows[226]["capability"]
    assert "complete naturalRemainder(x1)" in rows[226]["remaining_boundary"]
    assert "total-order-five axial derivatives" in rows[227]["capability"]
    assert "upstream contracts" in rows[227]["remaining_boundary"]
    assert "finite head" in rows[228]["capability"]
    assert "compact mean correction" in rows[228]["remaining_boundary"]
    assert "grad(Delta c)" in rows[229]["capability"]
    assert rows[229]["provenance"] is None
    assert "allJetsFlat_residual_of_stages" in rows[231]["capability"]
    assert "actual stage family" in rows[231]["remaining_boundary"]
    assert "hierarchy-owned" in rows[232]["capability"]
    assert "Issue-#1 leading fifth-mixed profile" in rows[232]["remaining_boundary"]
    assert "ProfileFourthMixedJet" in rows[233]["capability"]
    assert "explicit upstream contracts" in rows[233]["remaining_boundary"]


def test_cross_agent_source_truth_remains_fail_closed() -> None:
    pr226 = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2_axial_quadratic.json"
    )
    assert pr226["full_reconstruction"] is False
    assert pr226["paper_exact_velocity_available"] is False
    assert pr226["layer"]["status"] == "formal-structure"

    pr227 = _load(ROOT / "references" / "provenance_manifest_addendum_227.json")
    assert pr227["full_reconstruction"] is False
    assert pr227["paper_exact_velocity_available"] is False

    pr228 = _load(ROOT / "references" / "ACTUAL_SIGNED_MEAN_DEFECT_PROVENANCE.json")
    assert pr228["full_reconstruction"] is False
    assert pr228["paper_exact_velocity_available"] is False
    assert pr228["layer"]["truth_boundary"]["finite_head_mean_debt_materialized"] is False
    assert pr228["layer"]["truth_boundary"]["compact_mean_correction_available"] is False

    pr231 = _load(ROOT / "references" / "provenance_manifest_addendum_231.json")
    assert pr231["full_reconstruction"] is False
    assert pr231["paper_exact_velocity_available"] is False

    pr232 = _load(ROOT / "references" / "provenance_manifest_addendum_232.json")
    assert pr232["full_reconstruction"] is False
    assert pr232["paper_exact_velocity_available"] is False

    pr233 = _load(ROOT / "references" / "provenance_manifest_addendum_233.json")
    assert pr233["full_reconstruction"] is False
    assert pr233["paper_exact_velocity_available"] is False
    assert pr233["layers"][0]["status"] == "formal-structure"

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 4, 5, 6, 7):
        assert stages[stage_id]["status"] == "formal-structure"

    try:
        require_complete_reconstruction()
    except IncompleteReconstructionError:
        pass
    else:
        raise AssertionError("paper-exact completion gate must remain fail-closed")


def test_forbidden_inferences_and_agent_contract_stay_locked() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("slow2(x1)" in item and "x2" in item for item in forbidden)
    assert any("partial_eta^2 actualLowerSource/forcing" in item for item in forbidden)
    assert any("requested_cross_defect" in item and "mean debt" in item for item in forbidden)
    assert any("allJetsFlat_residual_of_stages" in item and "Section-9 stages" in item for item in forbidden)
    assert any("Section-10 cutoff" in item and "genuine forcing" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
