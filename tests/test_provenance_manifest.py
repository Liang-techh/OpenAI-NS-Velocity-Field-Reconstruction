import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "references" / "provenance_manifest.json"
ALLOWED = {"paper-exact", "formal-structure", "diagnostic-only", "toy", "pending"}


def load_manifest() -> dict:
    with MANIFEST.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def test_manifest_is_well_formed_and_source_pinned() -> None:
    data = load_manifest()
    assert data["schema_version"] == 1
    assert set(data["status_vocabulary"]) == ALLOWED

    commit = data["sources"]["official_lean"]["commit"]
    assert re.fullmatch(r"[0-9a-f]{40}", commit)
    assert data["sources"]["official_lean"]["repository"] == "openai/NavierStokesAndEuler"


def test_manifest_artifacts_exist_and_statuses_are_valid() -> None:
    data = load_manifest()
    ids = set()
    for layer in data["layers"]:
        assert layer["id"] not in ids
        ids.add(layer["id"])
        assert layer["status"] in ALLOWED
        assert layer["paper_locations"]
        for artifact in layer["artifacts"]:
            assert (ROOT / artifact).is_file(), artifact


def test_paper_exact_layers_have_no_recorded_blockers() -> None:
    data = load_manifest()
    for layer in data["layers"]:
        blockers = layer["missing_for_paper_exact"]
        if layer["status"] == "paper-exact":
            assert blockers == [], f"{layer['id']} is marked paper-exact but still has blockers"
        elif layer["status"] in {"formal-structure", "pending"}:
            assert blockers, f"{layer['id']} must record why it is not paper-exact"


def test_no_pending_layer_is_backed_by_claimed_runtime_artifacts() -> None:
    """A pending layer must not silently acquire implementation artifacts without review.

    When a stage starts being implemented, update its status to formal-structure (or stronger)
    in the same change that adds its manifest artifacts. This makes status drift review-visible.
    """
    data = load_manifest()
    for layer in data["layers"]:
        if layer["status"] == "pending":
            assert layer["artifacts"] == []


def test_manifest_tracks_latest_landed_formal_structure_without_promotion() -> None:
    data = load_manifest()
    layers = {layer["id"]: layer for layer in data["layers"]}

    assert data["full_reconstruction"] is False

    stage1 = layers["stage-1-leading-profile"]
    assert stage1["status"] == "formal-structure"
    assert "src/openai_ns_reconstruction/schedule_pressure.py" in stage1["artifacts"]
    assert "NavierStokes/OutgoingSchedule.lean" in stage1["official_lean_modules"]
    assert "NavierStokes/OutgoingTail.lean" in stage1["official_lean_modules"]
    assert any("S=32" in item and "flattenLength" in item for item in stage1["implemented_components"])
    assert any("finalAngular/clockWeight/axisPressure" in item for item in stage1["missing_for_paper_exact"])
    assert any("Classical.choose" in item for item in stage1["numerical_caveats"])

    stage2 = layers["stage-2-all-order-background"]
    assert stage2["status"] == "formal-structure"
    assert any("Eq. (5.6)" in item and "Omega_k/X" in item for item in stage2["implemented_components"])
    assert not any("construct the regular Omega" in item for item in stage2["missing_for_paper_exact"])

    stage3 = layers["stage-3-to-6-oscillatory-corrections"]
    assert stage3["status"] == "formal-structure"
    assert "src/openai_ns_reconstruction/phase_box_certificate.py" in stage3["artifacts"]
    assert any("rounded_normal_estimates" in item for item in stage3["implemented_components"])

    stage7 = layers["stage-7-final-localization"]
    assert stage7["status"] == "formal-structure"
    assert "src/openai_ns_reconstruction/time_localization.py" in stage7["artifacts"]
    assert "NavierStokes/TimeLocalization.lean" in stage7["paper_locations"]
    assert any("timeSwitch" in item for item in stage7["implemented_components"])
    assert any("smooth global force extension through t=1" in item for item in stage7["missing_for_paper_exact"])
    assert not any("paper-specific time localization/activation" in item for item in stage7["missing_for_paper_exact"])
