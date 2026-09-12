import json
from pathlib import Path

from openai_ns_reconstruction.provenance import (
    IncompleteReconstructionError,
    require_complete_reconstruction,
)
from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_237_integration.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_pr237_current_main_integration_stays_bounded() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "f3545b839ccd4c91dc207b4ca4c4020e8ec997e0"
    assert ledger["integrated_prs"] == [237]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    row = ledger["layers"][0]
    assert row["pr"] == 237
    assert row["stage"] == "stage-7-final-localization"
    assert row["status"] == "formal-structure"
    assert (ROOT / row["artifact"]).is_file()
    assert (ROOT / row["test"]).is_file()
    assert (ROOT / row["provenance"]).is_file()
    assert "off-axis spatial Laplacian" in row["capability"]
    assert "symmetry axis remains fail-closed" in row["remaining_boundary"]
    assert "genuine Navier-Stokes residual/forcing" in row["remaining_boundary"]


def test_pr237_does_not_open_paper_exact_gate() -> None:
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


def test_pr237_forbidden_promotions_remain_locked() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("axis regularity" in item for item in forbidden)
    assert any("genuine Navier-Stokes residual" in item for item in forbidden)
    assert any("finite-difference regression oracle" in item and "not a paper-exact proof" in item for item in forbidden)
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
