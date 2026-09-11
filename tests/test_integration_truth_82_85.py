import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_runtime_truth_tracks_prs_82_85_without_paper_exact_promotion() -> None:
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False
    assert stages[1]["status"] == "formal-structure"
    assert "axis_reference_pair" in stages[1]["implemented"]
    assert "zeroth parameter-jet radial coefficient functions" in stages[1]["implemented"]
    assert "full coefficient window" in stages[1]["implemented"]
    assert "complete compatible AxisSpace state" in stages[1]["implemented"]
    assert "Lift the landed actual referencePair radial coefficient functions" in stages[1]["remaining"]
    assert "coefficientOperators/naturalRemainder" in stages[1]["remaining"]

    assert stages[2]["status"] == "formal-structure"
    assert "lower-history bridge" in stages[2]["implemented"]
    assert "strict lower-order ProfileSecondJet history" in stages[2]["implemented"]
    assert "Z_(b-D) Z_b F_(n-1)" in stages[2]["implemented"]
    assert "Omega_(n-1)/X" in stages[2]["implemented"]
    assert "genuinely recursively solved/repaired orders < n" in stages[2]["remaining"]

    assert stages[4]["status"] == "formal-structure"
    assert "cylindrical_curl_jet" in stages[4]["implemented"]
    assert "B_theta/r" in stages[4]["implemented"]
    assert "|n|^-2(n cross a)" in stages[4]["implemented"]
    assert "actual localized paper coefficient" in stages[4]["implemented"]
    assert "cylindrical derivative jet accepted by cylindrical_curl_jet" in stages[4]["remaining"]

    assert stages[7]["status"] == "formal-structure"
    assert "section10_endpoint_ladder" in stages[7]["implemented"]
    assert "degrees 0..N" in stages[7]["implemented"]
    assert "same spatial window" in stages[7]["implemented"]
    assert "finite-order prerequisite interface" in stages[7]["implemented"]
    assert "close the infinite all-order family" in stages[7]["remaining"]

    limitations = " ".join(status["limitations"])
    assert "referencePair radial coefficient functions" in limitations
    assert "lower-history source adapter" in limitations
    assert "cylindrical_curl_jet" in limitations
    assert "finite section10_endpoint_ladder" in limitations


def test_machine_readable_addendum_maps_prs_82_85_and_remains_fail_closed() -> None:
    path = ROOT / "references" / "provenance_manifest_addendum_82_85.json"
    addendum = json.loads(path.read_text(encoding="utf-8"))

    assert addendum["kind"] == "provenance-manifest-addendum"
    assert addendum["base_manifest"] == "references/provenance_manifest.json"
    assert addendum["base_main_commit"] == "30da864f2f2ae1fc1d886c9c48e3a9763ceca658"
    assert addendum["integrated_prs"] == [82, 83, 84, 85]
    assert addendum["full_reconstruction"] is False
    assert addendum["paper_exact_velocity_available"] is False

    layers = {layer["id"]: layer for layer in addendum["layers"]}
    assert set(layers) == {
        "stage-1-leading-profile",
        "stage-2-all-order-background",
        "stage-3-to-6-oscillatory-corrections",
        "stage-7-final-localization",
    }
    assert all(layer["status"] == "formal-structure" for layer in layers.values())

    assert layers["stage-1-leading-profile"]["artifact"].endswith("axis_reference_pair.py")
    assert "naturalRemainder" in layers["stage-1-leading-profile"]["remaining_boundary"]
    assert layers["stage-2-all-order-background"]["artifact"].endswith("background_lower_history_source.py")
    assert "strict lower-order" in layers["stage-2-all-order-background"]["landed_capability"]
    assert layers["stage-3-to-6-oscillatory-corrections"]["artifact"].endswith("cylindrical_curl_jet.py")
    assert "B_theta/r" in layers["stage-3-to-6-oscillatory-corrections"]["landed_capability"]
    assert layers["stage-7-final-localization"]["artifact"].endswith("section10_endpoint_ladder.py")
    assert "infinite all-order" in layers["stage-7-final-localization"]["remaining_boundary"]


def test_readme_mentions_latest_integrated_boundaries_without_overclaiming() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "axis_reference_pair.py" in readme
    assert "background_lower_history_source.py" in readme
    assert "cylindrical_curl_jet.py" in readme
    assert "section10_endpoint_ladder.py" in readme
    assert "paper_exact_velocity_available: false" in readme
    assert "A successful demo or a green test suite is **not** a reconstruction" in readme
