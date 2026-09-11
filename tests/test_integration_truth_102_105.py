import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_prs_102_105_do_not_relax_runtime_completion_gate() -> None:
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False
    assert status["status"] == "partial-executable-reconstruction"
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[7]["status"] == "formal-structure"
    assert all(stage["status"] != "paper-exact" for stage in stages.values() if stage["id"] != 0)


def test_machine_readable_addendum_maps_prs_102_105_fail_closed() -> None:
    path = ROOT / "references" / "provenance_manifest_addendum_102_105.json"
    addendum = json.loads(path.read_text(encoding="utf-8"))

    assert addendum["kind"] == "provenance-manifest-addendum"
    assert addendum["base_manifest"] == "references/provenance_manifest.json"
    assert addendum["base_main_commit"] == "73ddba6a150bc932f8de7b533e74c80592cd2ef5"
    assert addendum["integrated_prs"] == [102, 103, 104, 105]
    assert addendum["full_reconstruction"] is False
    assert addendum["paper_exact_velocity_available"] is False

    layers = {layer["id"]: layer for layer in addendum["layers"]}
    assert set(layers) == {
        "stage-1-leading-profile-coefficient-product",
        "stage-2-all-order-background-third-mixed-repair",
        "stage-6-section9-flat-remainder-ladder",
        "stage-7-8-section9-endpoint-ladder-bridge",
    }
    assert all(layer["status"] == "formal-structure" for layer in layers.values())

    expected = {
        "stage-1-leading-profile-coefficient-product": (
            "src/openai_ns_reconstruction/axis_coefficient_product.py",
            "references/AXIS_COEFFICIENT_PRODUCT_PROVENANCE.md",
        ),
        "stage-2-all-order-background-third-mixed-repair": (
            "src/openai_ns_reconstruction/background_moment_repair_third_mixed_jets.py",
            "references/LEMMA_5_2_THIRD_MIXED_JET_PROVENANCE.md",
        ),
        "stage-6-section9-flat-remainder-ladder": (
            "src/openai_ns_reconstruction/section9_flat_remainder.py",
            "references/SECTION9_FLAT_REMAINDER_LADDER_PROVENANCE.md",
        ),
        "stage-7-8-section9-endpoint-ladder-bridge": (
            "src/openai_ns_reconstruction/section9_endpoint_ladder_bridge.py",
            "references/SECTION9_ENDPOINT_LADDER_BRIDGE_PROVENANCE.md",
        ),
    }
    for layer_id, (artifact, provenance) in expected.items():
        layer = layers[layer_id]
        assert layer["artifact"] == artifact
        assert layer["provenance"] == provenance
        assert (ROOT / artifact).is_file()
        assert (ROOT / provenance).is_file()
        assert layer["landed_capability"]
        assert layer["remaining_boundary"]

    assert "productFamily" in layers["stage-1-leading-profile-coefficient-product"]["landed_capability"]
    assert "naturalRemainder" in layers["stage-1-leading-profile-coefficient-product"]["remaining_boundary"]
    assert "U_XXeta" in layers["stage-2-all-order-background-third-mixed-repair"]["landed_capability"]
    assert "hierarchy-owned" in layers["stage-2-all-order-background-third-mixed-repair"]["remaining_boundary"]
    assert "N=0..Nmax" in layers["stage-6-section9-flat-remainder-ladder"]["landed_capability"]
    assert "all-orders flatness" in layers["stage-6-section9-flat-remainder-ladder"]["remaining_boundary"]
    assert "endpoint degrees 0..N" in layers["stage-7-8-section9-endpoint-ladder-bridge"]["landed_capability"]
    assert "finite q-flat ladder" in layers["stage-7-8-section9-endpoint-ladder-bridge"]["landed_capability"]
    assert "smooth compact forcing" in layers["stage-7-8-section9-endpoint-ladder-bridge"]["remaining_boundary"]


def test_one_off_pr_102_105_ledgers_and_readme_remain_fail_closed() -> None:
    one_off_addenda = (
        "references/provenance_manifest_addendum_axis_coefficient_product.json",
        "references/provenance_manifest_addendum_103.json",
        "references/provenance_manifest_addendum_104.json",
        "references/provenance_manifest_addendum_105.json",
    )
    for relative in one_off_addenda:
        one_off = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        assert one_off["full_reconstruction"] is False
        assert one_off["paper_exact_velocity_available"] is False

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "paper_exact_velocity_available: false" in readme
    assert "A successful demo or a green test suite is **not** a reconstruction" in readme
