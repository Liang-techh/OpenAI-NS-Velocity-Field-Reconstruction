import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_prs_97_100_do_not_relax_runtime_completion_gate() -> None:
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False
    assert status["status"] == "partial-executable-reconstruction"
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[7]["status"] == "formal-structure"
    assert all(stage["status"] != "paper-exact" for stage in stages.values() if stage["id"] != 0)


def test_machine_readable_addendum_maps_prs_97_100_fail_closed() -> None:
    path = ROOT / "references" / "provenance_manifest_addendum_97_100.json"
    addendum = json.loads(path.read_text(encoding="utf-8"))

    assert addendum["kind"] == "provenance-manifest-addendum"
    assert addendum["base_manifest"] == "references/provenance_manifest.json"
    assert addendum["base_main_commit"] == "5b851bbf187e9707e60cc79d59fb98a068c91150"
    assert addendum["integrated_prs"] == [97, 98, 99, 100]
    assert addendum["full_reconstruction"] is False
    assert addendum["paper_exact_velocity_available"] is False

    layers = {layer["id"]: layer for layer in addendum["layers"]}
    assert set(layers) == {
        "stage-1-leading-profile",
        "stage-2-all-order-background",
        "stage-6-section9-stage-certificate",
        "stage-7-8-endpoint-admission",
    }
    assert all(layer["status"] == "formal-structure" for layer in layers.values())

    expected = {
        "stage-1-leading-profile": (
            "src/openai_ns_reconstruction/axis_coefficient_reference_state.py",
            "references/AXIS_COEFFICIENT_REFERENCE_STATE_PROVENANCE.md",
        ),
        "stage-2-all-order-background": (
            "src/openai_ns_reconstruction/background_regular_flux_second_jets.py",
            "references/EQ_5_2_REGULAR_FLUX_SECOND_JET_PROVENANCE.md",
        ),
        "stage-6-section9-stage-certificate": (
            "src/openai_ns_reconstruction/section9_stage_certificate.py",
            "references/SECTION9_STAGE_CERTIFICATE_PROVENANCE.md",
        ),
        "stage-7-8-endpoint-admission": (
            "src/openai_ns_reconstruction/section9_endpoint_bridge.py",
            "references/SECTION9_ENDPOINT_BRIDGE_PROVENANCE.md",
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

    assert "coefficient-state view" in layers["stage-1-leading-profile"]["landed_capability"]
    assert "naturalRemainder" in layers["stage-1-leading-profile"]["remaining_boundary"]
    assert "beta_n=V_n/X" in layers["stage-2-all-order-background"]["landed_capability"]
    assert "U_XXeta" in layers["stage-2-all-order-background"]["remaining_boundary"]
    assert "Eq. (9.17)" in layers["stage-6-section9-stage-certificate"]["landed_capability"]
    assert "sum (9.21)" in layers["stage-6-section9-stage-certificate"]["remaining_boundary"]
    assert "order n+1" in layers["stage-7-8-endpoint-admission"]["landed_capability"]
    assert "t>=3/4" in layers["stage-7-8-endpoint-admission"]["landed_capability"]
    assert "uniform residual-derivative witness" in layers["stage-7-8-endpoint-admission"]["remaining_boundary"]
    assert "smooth force closure" in layers["stage-7-8-endpoint-admission"]["remaining_boundary"]


def test_one_off_pr_97_100_ledgers_and_readme_remain_fail_closed() -> None:
    one_off_addenda = (
        "references/provenance_manifest_addendum_axis_coefficient_reference_state.json",
        "references/provenance_manifest_addendum_98.json",
        "references/provenance_manifest_addendum_99.json",
        "references/provenance_manifest_addendum_100.json",
    )
    for relative in one_off_addenda:
        one_off = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        assert one_off["full_reconstruction"] is False
        assert one_off["paper_exact_velocity_available"] is False

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "paper_exact_velocity_available: false" in readme
    assert "A successful demo or a green test suite is **not** a reconstruction" in readme
