import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "references" / "provenance_manifest_addendum_112_115.json"


def _load(path: str | Path) -> dict:
    return json.loads(
        (ROOT / path).read_text(encoding="utf-8")
        if isinstance(path, str)
        else path.read_text(encoding="utf-8")
    )


def test_pr112_115_integration_ledger_is_fail_closed_and_resolves_to_landed_artifacts() -> None:
    ledger = _load(LEDGER_PATH)

    assert ledger["base_main_commit"] == "10889752d6bc7656d36f5578b26b61e27d2e5a11"
    assert ledger["integrated_prs"] == [112, 113, 114, 115]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False
    assert "remain formal-structure" in ledger["integration_boundary"]
    assert "paper_exact_velocity_available must remain false" in ledger["integration_boundary"]

    layers = {entry["pr"]: entry for entry in ledger["layers"]}
    assert set(layers) == {112, 113, 114, 115}

    for pr, entry in layers.items():
        assert entry["status"] == "formal-structure", pr
        for key in ("artifact", "test", "provenance", "one_off_manifest"):
            assert (ROOT / entry[key]).is_file(), (pr, key, entry[key])
        assert entry["remaining_boundary"], pr

    assert "complete coefficientOperators" in layers[112]["remaining_boundary"]
    assert "genuine naturalRemainder(x0)" in layers[112]["remaining_boundary"]
    assert "hierarchy moment/patch eta-jets remain inputs" in layers[113]["remaining_boundary"]
    assert "constructed Eq. (9.21) sums" in layers[114]["remaining_boundary"]
    assert "smooth extension through t=1" in layers[115]["remaining_boundary"]
    assert "smooth compact forcing" in layers[115]["remaining_boundary"]


def test_one_off_provenance_and_runtime_status_do_not_promote_pr112_115() -> None:
    ledger = _load(LEDGER_PATH)

    for entry in ledger["layers"]:
        one_off = _load(entry["one_off_manifest"])
        assert one_off["full_reconstruction"] is False, entry["pr"]
        assert one_off["paper_exact_velocity_available"] is False, entry["pr"]

        if "layer" in one_off:
            assert one_off["layer"]["status"] == "formal-structure", entry["pr"]
        else:
            assert one_off["layers"], entry["pr"]
            assert all(layer["status"] == "formal-structure" for layer in one_off["layers"]), entry["pr"]

    runtime = construction_status()
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert runtime["full_reconstruction"] is False
    assert runtime["paper_exact_velocity_available"] is False
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[6]["status"] != "paper-exact"
    assert stages[7]["status"] == "formal-structure"
