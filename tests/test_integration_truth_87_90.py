import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_prs_87_90_are_recorded_without_relaxing_completion_gates() -> None:
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False
    assert status["status"] == "partial-executable-reconstruction"
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[5]["status"] in {"pending", "formal-structure"}
    assert stages[7]["status"] == "formal-structure"

    # None of the newly landed artifacts is allowed to make the runtime gate
    # claim a completed paper instance merely by existing or passing tests.
    assert all(stage["status"] != "paper-exact" for stage in stages.values() if stage["id"] != 0)


def test_machine_readable_addendum_maps_prs_87_90_fail_closed() -> None:
    path = ROOT / "references" / "provenance_manifest_addendum_87_90.json"
    addendum = json.loads(path.read_text(encoding="utf-8"))

    assert addendum["kind"] == "provenance-manifest-addendum"
    assert addendum["base_manifest"] == "references/provenance_manifest.json"
    assert addendum["base_main_commit"] == "bdd7db01fe98e7eb216480e24384eadc7ccbf533"
    assert addendum["integrated_prs"] == [87, 88, 89, 90]
    assert addendum["full_reconstruction"] is False
    assert addendum["paper_exact_velocity_available"] is False

    layers = {layer["id"]: layer for layer in addendum["layers"]}
    assert set(layers) == {
        "stage-1-leading-profile",
        "stage-2-all-order-background",
        "stage-5-compact-mean-corrections",
        "stage-7-final-localization",
    }
    assert all(layer["status"] == "formal-structure" for layer in layers.values())

    expected = {
        "stage-1-leading-profile": (
            "src/openai_ns_reconstruction/axis_reference_angular_jets.py",
            "references/AXIS_REFERENCE_ANGULAR_JETS_PROVENANCE.md",
        ),
        "stage-2-all-order-background": (
            "src/openai_ns_reconstruction/background_lower_history_solver.py",
            "references/LOWER_HISTORY_FIRST_PICARD_PROVENANCE.md",
        ),
        "stage-5-compact-mean-corrections": (
            "src/openai_ns_reconstruction/mean_rank_update.py",
            "references/MEAN_RANK_UPDATE_PROVENANCE.md",
        ),
        "stage-7-final-localization": (
            "src/openai_ns_reconstruction/section10_endpoint_trace.py",
            "references/SECTION10_ENDPOINT_TRACE_PROVENANCE.md",
        ),
    }
    for layer_id, (artifact, provenance) in expected.items():
        layer = layers[layer_id]
        assert layer["artifact"] == artifact
        assert layer["provenance"] == provenance
        assert (ROOT / artifact).is_file()
        assert (ROOT / provenance).is_file()
        assert layer["remaining_boundary"]

    assert "arbitrary finite eta-derivative order" in layers["stage-1-leading-profile"]["landed_capability"]
    assert "axial reference coefficient" in layers["stage-1-leading-profile"]["remaining_boundary"]
    assert "xi to X=xi^2" in layers["stage-2-all-order-background"]["landed_capability"]
    assert "converged full Picard series" in layers["stage-2-all-order-background"]["remaining_boundary"]
    assert "five-row compact mean-repair" in layers["stage-5-compact-mean-corrections"]["landed_capability"]
    assert "actual mean debt" in layers["stage-5-compact-mean-corrections"]["remaining_boundary"]
    assert "endpoint enclosures" in layers["stage-7-final-localization"]["landed_capability"]
    assert "Actual derivative majorants" in layers["stage-7-final-localization"]["remaining_boundary"]


def test_new_provenance_notes_preserve_formal_structure_boundaries() -> None:
    notes = [
        "references/AXIS_REFERENCE_ANGULAR_JETS_PROVENANCE.md",
        "references/LOWER_HISTORY_FIRST_PICARD_PROVENANCE.md",
        "references/MEAN_RANK_UPDATE_PROVENANCE.md",
        "references/SECTION10_ENDPOINT_TRACE_PROVENANCE.md",
    ]
    for relative in notes:
        text = (ROOT / relative).read_text(encoding="utf-8")
        assert "formal-structure" in text
        assert "paper_exact_velocity_available" in text
        assert "false" in text.lower()

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "paper_exact_velocity_available: false" in readme
    assert "A successful demo or a green test suite is **not** a reconstruction" in readme
