import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_runtime_truth_surface_tracks_recent_landed_capabilities_without_promotion() -> None:
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert status["paper_exact_velocity_available"] is False
    assert stages[1]["status"] == stages[2]["status"] == stages[4]["status"] == stages[7]["status"] == "formal-structure"

    assert "actual SchedulePressure certificate is now wired" in stages[1]["implemented"]
    assert "positive factorial-majorant term" in stages[1]["implemented"]
    assert "binary64" in stages[1]["remaining"]

    assert "SlowExpansionResidual finite recurrence / truncation bridge" in stages[2]["implemented"]
    assert "first omitted slow factor" in stages[2]["implemented"]
    assert "independently proved identities" in stages[2]["remaining"]

    assert "actual-vs-reference 2x2 covariance perturbation certificate" in stages[4]["implemented"]
    assert "spectral ||H^-1||_2" in stages[4]["implemented"]
    assert "Eq. (7.28) analytic column-error bound" in stages[4]["remaining"]

    assert "sum_{j>N} 2^-j = 2^-N" in stages[7]["implemented"]
    assert "actual closed-past localized NS residual" in stages[7]["remaining"]


def test_manifest_maps_recent_modules_and_keeps_full_reconstruction_false() -> None:
    manifest = json.loads((ROOT / "references" / "provenance_manifest.json").read_text(encoding="utf-8"))
    layers = {layer["id"]: layer for layer in manifest["layers"]}

    assert manifest["full_reconstruction"] is False

    stage1 = layers["stage-1-leading-profile"]
    assert "src/openai_ns_reconstruction/stage1_scale_chain.py" in stage1["artifacts"]
    assert any("positive factorial-majorant term" in item and "binary64" in item
               for item in stage1["implemented_components"])
    assert any("not a lower bound on the actual natural resolvent" in item
               for item in stage1["numerical_caveats"])

    stage2 = layers["stage-2-all-order-background"]
    assert "src/openai_ns_reconstruction/background_truncation_residual.py" in stage2["artifacts"]
    assert any("recurrence_truncation" in item and "R_n" in item
               for item in stage2["implemented_components"])

    stage3 = layers["stage-3-to-6-oscillatory-corrections"]
    assert "src/openai_ns_reconstruction/stress_cone_perturbation.py" in stage3["artifacts"]
    assert any("determinant lower bound" in item and "spectral inverse bound" in item
               for item in stage3["implemented_components"])
    assert any("Eq. (7.28)" in item and "actual pulse-integrated covariance matrix H" in item
               for item in stage3["missing_for_paper_exact"])

    stage7 = layers["stage-7-final-localization"]
    assert any("DerivativeTailSumCertificate" in item and "2^-N" in item
               for item in stage7["implemented_components"])
    assert any("caller-supplied jets or analytic template bounds" in item and "2^-N" in item
               for item in stage7["numerical_caveats"])
