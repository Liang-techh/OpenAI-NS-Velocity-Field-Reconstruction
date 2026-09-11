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
    assert "landed wide Lambda/C selector" in stages[1]["implemented"]
    assert "symbolic exponent" in stages[1]["implemented"]
    assert "without binary64 down-conversion" in stages[1]["implemented"]
    assert "axis_fixed_point_picard" in stages[1]["implemented"]
    assert "s*remainderBound<=1" in stages[1]["implemented"]
    assert "s*remainderLip<=1/2" in stages[1]["implemented"]
    assert "s=1/(2Lambda)" in stages[1]["implemented"]
    assert "no coefficient-space fixed point has been materialized" in stages[1]["implemented"]
    assert "AxisCoefficientSpace/coefficientOperators/naturalRemainder backend" in stages[1]["remaining"]
    assert "pinned referencePair" in stages[1]["remaining"]
    assert "materialize the coefficient-space fixed-point phi/u" in stages[1]["remaining"]
    assert "former binary64 Lambda/C representation blocker is closed" in stages[1]["remaining"]
    assert "not a lower bound" in stages[1]["remaining"]

    assert "SlowExpansionResidual finite recurrence / truncation bridge" in stages[2]["implemented"]
    assert "finite eta-jet" in stages[2]["implemented"]
    assert "p=e_*f" in stages[2]["implemented"]
    assert "function-level compact-repair adapter" in stages[2]["implemented"]
    assert "paper_exact=False" in stages[2]["implemented"]
    assert "arbitrary-precision power-of-two Python integer witness" in stages[2]["implemented"]
    assert "reciprocal_support_log_edge" in stages[2]["implemented"]
    assert "first omitted slow factor" in stages[2]["implemented"]
    assert "PositiveAxisSystem adapter" in stages[2]["implemented"]
    assert "BaseJet(phi,axial,beta)" in stages[2]["implemented"]
    assert "SourceJet(angular,axial,pressureProduct,omegaQuotient)" in stages[2]["implemented"]
    assert "lambda_n=2nh" in stages[2]["implemented"]
    assert "materialized leading and lower-order BaseJet/SourceJet" in stages[2]["remaining"]
    assert "genuine moment/p=e_*f jets" in stages[2]["remaining"]
    assert "independently proved" in stages[2]["remaining"]

    assert "actual-vs-reference 2x2 covariance" in stages[4]["implemented"]
    assert "spectral ||H^-1||_2" in stages[4]["implemented"]
    assert "pulse-covariance error-budget adapter" in stages[4]["implemented"]
    assert "|e_sigma|<=E_ratio+u_*sqrt(1+c_*^2)M1" in stages[4]["implemented"]
    assert "coefficient-level curl-realization algebra" in stages[4]["implemented"]
    assert "Exact tangency n.a=0" in stages[4]["implemented"]
    assert "primary amplitude ODE algebra" in stages[4]["implemented"]
    assert "x=p+q, y=h(p-q)" in stages[4]["implemented"]
    assert "finite-interval PrimaryODE adapter" in stages[4]["implemented"]
    assert "Volterra integral defect" in stages[4]["implemented"]
    assert "closed-form constant-coefficient oracle" in stages[4]["implemented"]
    assert "not the paper's Proposition 5.5/pulse instance" in stages[4]["implemented"]
    assert "actual pulse-integrated covariance columns of Eq. (7.27)" in stages[4]["remaining"]
    assert "Gaussian normalized first-moment" in stages[4]["remaining"]
    assert "finite-interval PrimaryODE Volterra solution adapter" in stages[4]["remaining"]
    assert "rigorous error certificate" in stages[4]["remaining"]
    assert "actual localized coefficient" in stages[4]["remaining"]

    assert "sum_{j>N} 2^-j = 2^-N" in stages[7]["implemented"]
    assert "close_left_open_past" in stages[7]["implemented"]
    assert "rejects t>T" in stages[7]["implemented"]
    assert "endpoint-limit majorant" in stages[7]["implemented"]
    assert "0<=alpha<1" in stages[7]["implemented"]
    assert "endpoint-localization transfer" in stages[7]["implemented"]
    assert "requires T=1" in stages[7]["implemented"]
    assert "start time at or after 3/4" in stages[7]["implemented"]
    assert "timeSwitch'=0" in stages[7]["implemented"]
    assert "Transition-collar or non-T=1 majorants fail closed" in stages[7]["implemented"]
    assert "pointwise endpoint-force support certificate" in stages[7]["implemented"]
    assert "requiring exact zeros with no tolerance" in stages[7]["implemented"]
    assert "actual closed-past localized NS residual" in stages[7]["remaining"]
    assert "official t>=3/4 plateau" in stages[7]["remaining"]
    assert "fail-closed T=1 transfer" in stages[7]["remaining"]
    assert "degree-zero endpoint value" in stages[7]["remaining"]
    assert "every genuine endpoint normal coefficient vanish outside" in stages[7]["remaining"]


def test_manifest_maps_recent_modules_and_keeps_full_reconstruction_false() -> None:
    manifest = json.loads((ROOT / "references" / "provenance_manifest.json").read_text(encoding="utf-8"))
    layers = {layer["id"]: layer for layer in manifest["layers"]}

    assert manifest["full_reconstruction"] is False

    stage1 = layers["stage-1-leading-profile"]
    assert "src/openai_ns_reconstruction/axis_componentwise_input_bounds.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/axis_remainder_wide_bounds.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/natural_scale_selection_wide.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/stage1_scale_chain.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/stage1_scale_chain_wide.py" in stage1["artifacts"]
    assert "src/openai_ns_reconstruction/axis_fixed_point_picard.py" in stage1["artifacts"]
    assert any("chi-specific" in item and "AxisResolvent" in item
               for item in stage1["implemented_components"])
    assert any("factorial resolvent majorant is now representable" in item
               for item in stage1["implemented_components"])
    assert any("96-digit Decimal" in item and "sys.float_info.max" in item
               for item in stage1["implemented_components"])
    assert any("wide theorem-shaped Lambda/C continuation" in item and "SymbolicExponentialThreshold" in item
               for item in stage1["implemented_components"])
    assert any("AxisContraction scalar gate" in item and "s*remainderBound<=1" in item
               and "naturalRemainder" in item for item in stage1["implemented_components"])
    assert any("AxisCoefficientSpace" in item and "pinned referencePair" in item
               for item in stage1["missing_for_paper_exact"])
    assert any("materialize the coefficient-space fixed-point fields phi/u" in item
               for item in stage1["missing_for_paper_exact"])
    assert not any("through the pinned Lambda/C scale selection without binary64 down-conversion" in item
                   for item in stage1["missing_for_paper_exact"])
    assert any("axis_remainder_wide_bounds.py" in item and "not interval arithmetic" in item
               for item in stage1["numerical_caveats"])
    assert any("natural_scale_selection_wide.py" in item and "does not materialize" in item
               for item in stage1["numerical_caveats"])
    assert any("axis_fixed_point_picard.py" in item and "plumbing-only" in item
               for item in stage1["numerical_caveats"])

    stage2 = layers["stage-2-all-order-background"]
    assert "src/openai_ns_reconstruction/background_positive_axis.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_moment_repair_jets.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_moment_repair_profile.py" in stage2["artifacts"]
    assert "src/openai_ns_reconstruction/background_truncation_residual.py" in stage2["artifacts"]
    assert any("PositiveAxisSystem coefficient adapter" in item and "BaseJet" in item and "SourceJet" in item
               for item in stage2["implemented_components"])
    assert any("finite eta-jet" in item and "p=e_*f" in item
               for item in stage2["implemented_components"])
    assert any("function-level Lemma 5.2 compact repair" in item and "paper_exact" in item
               for item in stage2["implemented_components"])
    assert any("arbitrary-precision power-of-two integer witness" in item and "reciprocal_support_log_edge" in item
               for item in stage2["implemented_components"])
    assert any("recurrence_truncation" in item and "R_n" in item
               for item in stage2["implemented_components"])
    assert any("materialized leading BaseJet" in item and "lower-order SourceJet" in item
               for item in stage2["missing_for_paper_exact"])
    assert any("background_positive_axis.py" in item and "not the materialized paper hierarchy" in item
               for item in stage2["numerical_caveats"])
    assert any("background_moment_repair_jets.py" in item and "does not derive" in item
               for item in stage2["numerical_caveats"])
    assert any("background_moment_repair_profile.py" in item and "paper_exact remains False" in item
               for item in stage2["numerical_caveats"])
    assert any("arbitrary-precision integer/log-edge path" in item and "machine-representation obstruction" in item
               for item in stage2["numerical_caveats"])

    stage3 = layers["stage-3-to-6-oscillatory-corrections"]
    assert "src/openai_ns_reconstruction/stress_cone_perturbation.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/pulse_covariance_budget.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/curl_realization_algebra.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/primary_amplitude_ode.py" in stage3["artifacts"]
    assert "src/openai_ns_reconstruction/primary_amplitude_solution.py" in stage3["artifacts"]
    assert "NavierStokes/CurlClassBounds.lean" in stage3["official_lean_modules"]
    assert "NavierStokes/LocalizedCurlRealization.lean" in stage3["official_lean_modules"]
    assert "NavierStokes/MovingFrameODE.lean" in stage3["official_lean_modules"]
    assert "NavierStokes/GrowingMode.lean" in stage3["official_lean_modules"]
    assert "NavierStokes/PrimaryODE.lean" in stage3["official_lean_modules"]
    assert any("determinant lower bound" in item and "spectral inverse bound" in item
               for item in stage3["implemented_components"])
    assert any("pulse-covariance error-budget" in item and "C/sqrt(S_*)" in item
               for item in stage3["implemented_components"])
    assert any("Formula (30)" in item and "inverseCarrier=i/K" in item
               for item in stage3["implemented_components"])
    assert any("primary amplitude ODE pointwise algebra" in item and "modalOperator" in item
               for item in stage3["implemented_components"])
    assert any("finite-interval PrimaryODE execution adapter" in item and "Volterra defect" in item
               for item in stage3["implemented_components"])
    assert any("actual pulse-integrated covariance matrix H" in item and "Eq. (7.28)" in item
               for item in stage3["missing_for_paper_exact"])
    assert any("finite-interval PrimaryODE Volterra solution adapter" in item and "actual amplitude" in item
               for item in stage3["missing_for_paper_exact"])
    assert any("pulse_covariance_budget.py" in item and "does not construct the actual pulse" in item
               for item in stage3["numerical_caveats"])
    assert any("curl_realization_algebra.py" in item and "caller-supplied coefficient_curl" in item
               for item in stage3["numerical_caveats"])
    assert any("primary_amplitude_ode.py" in item and "does not construct the finite-interval Volterra solution" in item
               for item in stage3["numerical_caveats"])
    assert any("primary_amplitude_solution.py" in item and "closed-form regression oracle" in item
               and "does not instantiate" in item for item in stage3["numerical_caveats"])

    stage7 = layers["stage-7-final-localization"]
    assert "src/openai_ns_reconstruction/endpoint_limit_majorant.py" in stage7["artifacts"]
    assert "src/openai_ns_reconstruction/endpoint_support.py" in stage7["artifacts"]
    assert any("close_left_open_past" in item and "rejects t>T" in item
               for item in stage7["implemented_components"])
    assert any("DerivativeTailSumCertificate" in item and "2^-N" in item
               for item in stage7["implemented_components"])
    assert any("endpoint residual-jet limit majorant" in item and "alpha<1" in item
               for item in stage7["implemented_components"])
    assert any("endpoint-localization transfer" in item and "T=1" in item and "start>=3/4" in item
               for item in stage7["implemented_components"])
    assert any("pointwise endpoint-force support implication" in item and "no tolerance" in item
               for item in stage7["implemented_components"])
    assert any("endpoint-limit majorant" in item and "T=1 localization transfer" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("degree-zero endpoint value used by close_left_open_past" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("actual localized past residual" in item and "global compact-support conclusion" in item
               for item in stage7["missing_for_paper_exact"])
    assert any("endpoint_borel.close_left_open_past" in item and "does not prove continuity" in item
               for item in stage7["numerical_caveats"])
    assert any("endpoint_limit_majorant.py" in item and "finite samples" in item
               for item in stage7["numerical_caveats"])
    assert any("section10_endpoint_localization_transfer" in item and "does not derive" in item
               for item in stage7["numerical_caveats"])
    assert any("endpoint_support.py" in item and "finitely many point queries" in item
               for item in stage7["numerical_caveats"])
    assert any("caller-supplied jets or analytic template bounds" in item and "2^-N" in item
               for item in stage7["numerical_caveats"])
