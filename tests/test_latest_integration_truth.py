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
    assert "componentwise refinement" in stages[1]["implemented"]
    assert "chi-specific ||chi||<=12B_chi" in stages[1]["implemented"]
    assert "factorial resolvent majorant is representable" in stages[1]["implemented"]
    assert "96-digit Decimal" in stages[1]["implemented"]
    assert "exceeds sys.float_info.max" in stages[1]["implemented"]
    assert "remainderBound/remainderLip" in stages[1]["remaining"]
    assert "Lambda/C" in stages[1]["remaining"]
    assert "not a lower bound" in stages[1]["remaining"]

    assert "SlowExpansionResidual finite recurrence / truncation bridge" in stages[2]["implemented"]
    assert "finite eta-jet" in stages[2]["implemented"]
    assert "p=e_*f" in stages[2]["implemented"]
    assert "function-level compact-repair adapter" in stages[2]["implemented"]
    assert "paper_exact=False" in stages[2]["implemented"]
    assert "first omitted slow factor" in stages[2]["implemented"]
    assert "genuine moment/p=e_*f jets" in stages[2]["remaining"]
    assert "independently proved" in stages[2]["remaining"]

    assert "actual-vs-reference 2x2 covariance" in stages[4]["implemented"]
    assert "spectral ||H^-1||_2" in stages[4]["implemented"]
    assert "pulse-covariance error-budget adapter" in stages[4]["implemented"]
    assert "|e_sigma|<=E_ratio+u_*sqrt(1+c_*^2)M1" in stages[4]["implemented"]
    assert "coefficient-level curl-realization algebra" in stages[4]["implemented"]
    assert "Exact tangency n.a=0" in stages[4]["implemented"]
    assert "actual pulse-integrated covariance columns of Eq. (7.27)" in stages[4]["remaining"]
    assert "Gaussian normalized first-moment" in stages[4]["remaining"]
    assert "actual localized coefficient" in stages[4]["remaining"]

    assert "sum_{j>N} 2^-j = 2^-N" in stages[7]["implemented"]
    assert "close_left_open_past" in stages[7]["implemented"]
    assert "rejects t>T" in stages[7]["implemented"]
    assert "endpoint-limit majorant" in stages[7]["implemented"]
    assert "0<=alpha<1" in stages[7]["implemented"]
    assert "actual closed-past localized NS residual" in stages[7]["remaining"]
    assert "degree-zero endpoint value" in stages[7]["remaining"]
    assert "discharge the landed endpoint-limit majorant" in stages[7]["remaining"]


def test_manifest_maps_recent_modules_and_keeps_full_reconstruction_false() -> None:
    manifest = json.loads((ROOT / "references" / "provenance_manifest.json").read_text(encoding="utf-8"))
    layers = {layer["id"]: layer for layer in manifest["layers"]}

    assert manifest["full_reconstruction"] is False

    stage1 = layers["stage-1-leading-profile"]
    assert "src/openai_ns_reconstruction/axis_componentwise_input_bounds.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/axis_remainder_wide_bounds.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/stage1_scale_chain.py" in stage1["artifacts"]
    assert any("chi-specific" in item and "AxisResolvent" in item
               for item in stage1["implemented_components"])
    assert any("factorial resolvent majorant is now representable" in item
               for item in stage1["implemented_components"])
    assert any("96-digit Decimal" in item and "sys.float_info.max" in item
               for item in stage1["implemented_components"])
    assert any("remainderBound/remainderLip" in item and "not a lower bound" in item
               for item in stage1["missing_for_paper_exact"])
    assert any("axis_remainder_wide_bounds.py" in item and "not interval arithmetic" in item
               for item in stage1["numerical_caveats"])

    stage2 = layers["stage-2-all-order-background"]
    assert "src/openai_ns_reconstruction/background_moment_repair_jets.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_moment_repair_profile.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_truncation_residual.py" in stage2["artifacts"]
    assert any("finite eta-jet" in item and "p=e_*f" in item
               for item in stage2["implemented_components"])
    assert any("function-level Lemma 5.2 compact repair" in item and "paper_exact" in item
               for item in stage2["implemented_components"])
    assert any("recurrence_truncation" in item and "R_n" in item
               for item in stage2["implemented_components"])
    assert any("background_moment_repair_jets.py" in item and "does not derive" in item
               for item in stage2["numerical_caveats"])
    assert any("background_moment_repair_profile.py" in item and "paper_exact remains False" in item
               for item in stage2["numerical_caveats"])

    stage3 = layers["stage-3-to-6-oscillatory-corrections"]
    assert "src/openai_ns_reconstruction/stress_cone_perturbation.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/pulse_covariance_budget.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/curl_realization_algebra.py" in stage3["artifacts"]
    assert "NavierStokes/CurlClassBounds.lean" in stage3["official_lean_modules"]
    assert "NavierStokes/LocalizedCurlRealization.lean" in stage3["official_lean_modules"]
    assert any("determinant lower bound" in item and "spectral inverse bound" in item
               for item in stage3["implemented_components"])
    assert any("pulse-covariance error-budget" in item and "C/sqrt(S_*)" in item
               for item in stage3["implemented_components"])
    assert any("Formula (30)" in item and "inverseCarrier=i/K" in item
               for item in stage3["implemented_components"])
    assert any("actual pulse-integrated covariance matrix H" in item and "Eq. (7.28)" in item
               for item in stage3["missing_for_paper_exact"])
    assert any("actual localized coefficient" in item and "LocalizedCurlRealization" in item
               for item in stage3["missing_for_paper_exact"])
    assert any("pulse_covariance_budget.py" in item and "does not construct the actual pulse" in item
               for item in stage3["numerical_caveats"])
    assert any("curl_realization_algebra.py" in item and "caller-supplied coefficient_curl" in item
               for item in stage3["numerical_caveats"])

    stage7 = layers["stage-7-final-localization"]
    assert "src/openai_ns_reconstruction/endpoint_limit_majorant.py" in stage7["artifacts"]
    assert any("close_left_open_past" in item and "rejects t>T" in item
               for item in stage7["implemented_components"])
    assert any("DerivativeTailSumCertificate" in item and "2^-N" in item
               for item in stage7["implemented_components"])
    assert any("endpoint residual-jet limit majorant" in item and "alpha<1" in item
               for item in stage7["implemented_components"])
    assert any("endpoint-limit majorant" in item and "locally uniform limits" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("degree-zero endpoint value used by close_left_open_past" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("endpoint_borel.close_left_open_past" in item and "does not prove continuity" in item
               for item in stage7["numerical_caveats"])
    assert any("endpoint_limit_majorant.py" in item and "finite samples" in item
               for item in stage7["numerical_caveats"])
    assert any("caller-supplied jets or analytic template bounds" in item and "2^-N" in item
               for item in stage7["numerical_caveats"])
