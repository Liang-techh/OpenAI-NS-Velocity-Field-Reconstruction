import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "references" / "provenance_manifest_addendum_107_110.json"


def _load(path: str | Path) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8") if isinstance(path, str) else path.read_text(encoding="utf-8"))


def test_pr107_110_integration_ledger_is_fail_closed_and_resolves_to_landed_artifacts() -> None:
    ledger = _load(LEDGER_PATH)

    assert ledger["base_main_commit"] == "87e12b1cd5e1d973c4ada03fc9125cf7e7f8a74d"
    assert ledger["integrated_prs"] == [107, 108, 109, 110]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False
    assert "remain formal-structure" in ledger["integration_boundary"]
    assert "paper_exact_velocity_available must remain false" in ledger["integration_boundary"]

    layers = {entry["pr"]: entry for entry in ledger["layers"]}
    assert set(layers) == {107, 108, 109, 110}

    for pr, entry in layers.items():
        assert entry["status"] == "formal-structure", pr
        for key in ("artifact", "test", "provenance", "one_off_manifest"):
            assert (ROOT / entry[key]).is_file(), (pr, key, entry[key])
        assert entry["remaining_boundary"], pr

    assert "complete coefficientOperators" in layers[107]["remaining_boundary"]
    assert "genuine naturalRemainder(x0)" in layers[107]["remaining_boundary"]
    assert "hierarchy-owned moment/patch eta-jets remain inputs" in layers[108]["remaining_boundary"]
    assert "constructed Eq. (9.21) sums" in layers[109]["remaining_boundary"]
    assert "rejects t_upper=1" in layers[110]["remaining_boundary"]
    assert "uniform endpoint derivative majorants" in layers[110]["remaining_boundary"]


def test_one_off_provenance_and_runtime_status_do_not_promote_pr107_110() -> None:
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
    assert runtime["paper_exact_velocity_available"] is False
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[5]["status"] != "paper-exact"
    assert stages[6]["status"] != "paper-exact"
    assert stages[7]["status"] == "formal-structure"
