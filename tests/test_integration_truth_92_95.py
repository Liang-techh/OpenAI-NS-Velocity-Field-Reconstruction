import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_prs_92_95_do_not_relax_runtime_completion_gate() -> None:
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False
    assert status["status"] == "partial-executable-reconstruction"
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[7]["status"] == "formal-structure"
    assert all(stage["status"] != "paper-exact" for stage in stages.values() if stage["id"] != 0)


def test_machine_readable_addendum_maps_prs_92_95_fail_closed() -> None:
    path = ROOT / "references" / "provenance_manifest_addendum_92_95.json"
    addendum = json.loads(path.read_text(encoding="utf-8"))

    assert addendum["kind"] == "provenance-manifest-addendum"
    assert addendum["base_manifest"] == "references/provenance_manifest.json"
    assert addendum["base_main_commit"] == "bdba95a806f3ccb243956ff5873812b59a744011"
    assert addendum["integrated_prs"] == [92, 93, 94, 95]
    assert addendum["full_reconstruction"] is False
    assert addendum["paper_exact_velocity_available"] is False

    layers = {layer["id"]: layer for layer in addendum["layers"]}
    assert set(layers) == {
        "stage-1-leading-profile",
        "stage-2-all-order-background",
        "stage-6-section9-residual-improvement",
        "stage-7-final-localization",
    }
    assert all(layer["status"] == "formal-structure" for layer in layers.values())

    expected = {
        "stage-1-leading-profile": (
            "src/openai_ns_reconstruction/axis_reference_axial_jets.py",
            "references/AXIS_REFERENCE_AXIAL_JETS_PROVENANCE.md",
        ),
        "stage-2-all-order-background": (
            "src/openai_ns_reconstruction/background_moment_repair_second_jets.py",
            "references/LEMMA_5_2_SECOND_JET_PROVENANCE.md",
        ),
        "stage-6-section9-residual-improvement": (
            "src/openai_ns_reconstruction/section9_residual_decay.py",
            "references/SECTION9_RESIDUAL_DECAY_PROVENANCE.md",
        ),
        "stage-7-final-localization": (
            "src/openai_ns_reconstruction/spatial_localization.py",
            "references/SECTION10_FIXED_CUTOFF_BINDING_PROVENANCE.md",
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

    assert "axial reference coefficient" in layers["stage-1-leading-profile"]["landed_capability"]
    assert "naturalRemainder" in layers["stage-1-leading-profile"]["remaining_boundary"]
    assert "second (X,eta) jets" in layers["stage-2-all-order-background"]["landed_capability"]
    assert "beta_n=V_n/X" in layers["stage-2-all-order-background"]["remaining_boundary"]
    assert "rational arithmetic" in layers["stage-6-section9-residual-improvement"]["landed_capability"]
    assert "Proposition 9.9 residual flatness" in layers["stage-6-section9-residual-improvement"]["remaining_boundary"]
    assert "section10_localized_field" in layers["stage-7-final-localization"]["landed_capability"]
    assert "smooth compact forcing" in layers["stage-7-final-localization"]["remaining_boundary"]


def test_landed_provenance_and_existing_one_off_addenda_remain_fail_closed() -> None:
    notes = [
        "references/AXIS_REFERENCE_AXIAL_JETS_PROVENANCE.md",
        "references/LEMMA_5_2_SECOND_JET_PROVENANCE.md",
        "references/SECTION9_RESIDUAL_DECAY_PROVENANCE.md",
        "references/SECTION10_FIXED_CUTOFF_BINDING_PROVENANCE.md",
    ]
    for relative in notes:
        text = (ROOT / relative).read_text(encoding="utf-8")
        assert "formal-structure" in text
        assert "paper_exact_velocity_available" in text
        assert "false" in text.lower()

    for relative in (
        "references/provenance_manifest_addendum_axis_reference_axial_jets.json",
        "references/provenance_manifest_addendum_section9_residual_decay.json",
    ):
        one_off = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        assert one_off["full_reconstruction"] is False
        assert one_off["paper_exact_velocity_available"] is False

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "paper_exact_velocity_available: false" in readme
    assert "A successful demo or a green test suite is **not** a reconstruction" in readme
