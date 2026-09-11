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
    assert "src/openai_ns_reconstruction/outgoing_tail.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/schedule_axis_pressure.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/schedule_axis_margin.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/natural_scale_selection.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/axis_remainder_bounds.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/axis_analytic_input_bounds.py" in stage1["artifacts"]
    assert "NavierStokes/AxisContraction.lean" in stage1["official_lean_modules"]
    assert "NavierStokes/AxisEvaluation.lean" in stage1["official_lean_modules"]
    assert "NavierStokes/AnalyticCoefficientBounds.lean" in stage1["official_lean_modules"]
    assert "NavierStokes/AxisResolvent.lean" in stage1["official_lean_modules"]
    assert any("S=32" in item and "flattenLength" in item for item in stage1["implemented_components"])
    assert any("finalAngular" in item and "clockWeight" in item for item in stage1["implemented_components"])
    assert any("actual SchedulePressure.axisPressure" in item for item in stage1["implemented_components"])
    assert any("low-|Z|" in item and "sigma=sqrt(m)/20" in item for item in stage1["implemented_components"])
    assert any("Lambda=max(1+B+L,1+(14000/9)B)" in item and "C=exp" in item
               for item in stage1["implemented_components"])
    assert any("remainderBound/remainderLip" in item for item in stage1["implemented_components"])
    assert any("radiusLoss(1/2)=12" in item and "M=12" in item and "K<=30720B" in item
               for item in stage1["implemented_components"])
    assert any("analytic neighborhood radius rho" in item and "complex-field sup bound B" in item
               and "realPartSup" in item for item in stage1["missing_for_paper_exact"])
    assert not any("naturalResolvent factorial-series norm bound" in item
                   for item in stage1["missing_for_paper_exact"])
    assert any("Classical.choose" in item for item in stage1["numerical_caveats"])
    assert any("binary64" in item and "interval-arithmetic" in item for item in stage1["numerical_caveats"])
    assert any("axis_analytic_input_bounds.py" in item and "does not manufacture" in item
               for item in stage1["numerical_caveats"])

    stage2 = layers["stage-2-all-order-background"]
    assert stage2["status"] == "formal-structure"
    assert "src/openai_ns_reconstruction/background_inner_solver.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_moment_repair.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_extension.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_cutoff_schedule.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_cutoff_support.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_tail_order.py" in stage2["artifacts"]
    assert "Eq. (5.7)" in stage2["paper_locations"]
    assert any("Eq. (5.6)" in item and "Omega_k/X" in item for item in stage2["implemented_components"])
    assert any("Eq. (5.7)" in item and "Picard" in item for item in stage2["implemented_components"])
    assert any("Eq. (5.15) forward reconstruction" in item for item in stage2["implemented_components"])
    assert any("SlowBorelBase/DiagonalScale" in item and "C[j,m]" in item
               for item in stage2["implemented_components"])
    assert any("2^-J" in item and "h(J+1)-m" in item and "h(J+1)+b-2M" in item
               for item in stage2["implemented_components"])
    assert any("A0/A1/f_n" in item for item in stage2["missing_for_paper_exact"])
    assert any("true eta-dependent repaired coefficient hierarchy" in item
               for item in stage2["missing_for_paper_exact"])
    assert any("true uniform compactness bounds C[j,m]" in item
               for item in stage2["missing_for_paper_exact"])
    assert any("all-jets-flat residual decay" in item for item in stage2["missing_for_paper_exact"])
    assert any("background_tail_order.py" in item and "does not prove" in item
               for item in stage2["numerical_caveats"])

    stage3 = layers["stage-3-to-6-oscillatory-corrections"]
    assert stage3["status"] == "formal-structure"
    assert "src/openai_ns_reconstruction/slow_labels.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/phase_box_certificate.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/slot_geometry.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/slow_support_adjacency.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/phase_uniform_bounds.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/phase_frame_bounds.py" in stage3["artifacts"]
    assert "NavierStokes/BasePhaseGeometry.lean" in stage3["official_lean_modules"]
    assert any("Eq. (6.8)" in item for item in stage3["implemented_components"])
    assert any("2250-color" in item and "common r0" in item for item in stage3["implemented_components"])
    assert any("physical slow-support to SlotColoring.Adj bridge" in item
               for item in stage3["implemented_components"])
    assert any("rounded_normal_estimates" in item for item in stage3["implemented_components"])
    assert any("phaseConstant(M)=normalConstant(frequencyBound(M))" in item
               and "16 G^2(1+G)E" in item and "4 M(2A+5)delta" in item
               for item in stage3["implemented_components"])
    assert any("paper-exact Proposition 5.5 background" in item
               for item in stage3["missing_for_paper_exact"])
    assert any("LocalBaseBounds C1/C2" in item and "normal-closeness" in item
               for item in stage3["missing_for_paper_exact"])
    assert any("phase_frame_bounds.py" in item and "does not certify" in item
               for item in stage3["numerical_caveats"])

    stage7 = layers["stage-7-final-localization"]
    assert stage7["status"] == "formal-structure"
    assert "src/openai_ns_reconstruction/time_localization.py" in stage7["artifacts"]
    assert "src/openai_ns_reconstruction/past_extension.py" in stage7["artifacts"]
    assert "src/openai_ns_reconstruction/endpoint_borel.py" in stage7["artifacts"]
    assert "src/openai_ns_reconstruction/endpoint_scale_schedule.py" in stage7["artifacts"]
    assert "src/openai_ns_reconstruction/traced_residual.py" in stage7["artifacts"]
    assert "src/openai_ns_reconstruction/section10_support_energy.py" in stage7["artifacts"]
    assert "NavierStokes/PastExtension.lean" in stage7["paper_locations"]
    assert any("timeSwitch" in item for item in stage7["implemented_components"])
    assert any("PastExtension closed-past branch" in item for item in stage7["implemented_components"])
    assert any("doublingEnvelope" in item and "right-extension" in item for item in stage7["implemented_components"])
    assert any("boundSum/localScale" in item and "2^-j" in item for item in stage7["implemented_components"])
    assert any("pi/32" in item and "E(t)<=pi M^2/64" in item
               for item in stage7["implemented_components"])
    assert any("closed-past localized NS residual" in item and "full spacetime derivative family" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("analytic compact-template derivative bounds" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("smooth global force extension through t=1" in item for item in stage7["missing_for_paper_exact"])
    assert any("uniform bounded kinetic-energy" in item for item in stage7["missing_for_paper_exact"])
    assert any("tautological" in item for item in stage7["numerical_caveats"])
    assert any("section10_support_energy.py" in item and "fixed-time implication" in item
               for item in stage7["numerical_caveats"])
