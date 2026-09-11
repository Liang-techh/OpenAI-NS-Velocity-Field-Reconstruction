import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_122_125.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_landed_artifacts_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "fefcaf368d9c36b2b4bd48842c3b3d39e9f7d33d"
    assert ledger["integrated_prs"] == [122, 123, 124, 125]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {122, 123, 124, 125}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["remaining_boundary"].strip()

    assert "naturalRemainder(x0)" in rows[122]["remaining_boundary"]
    assert "partial_eta f_n remains an explicit analytic input" in rows[123]["remaining_boundary"]
    assert "sufficient gate" in rows[124]["remaining_boundary"]
    assert "No actual Section 9 residual jets" in rows[125]["remaining_boundary"]


def test_source_addenda_and_runtime_remain_fail_closed() -> None:
    addenda = (
        "provenance_manifest_addendum_axis_coefficient_multiply_y.json",
        "provenance_manifest_addendum_123.json",
        "provenance_manifest_addendum_phase_large_band_scale.json",
        "provenance_manifest_addendum_section10_endpoint_trace_consistency.json",
    )
    for name in addenda:
        data = _load(ROOT / "references" / name)
        assert data["full_reconstruction"] is False
        assert data["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[1]["status"] == "formal-structure"
    assert stages[2]["status"] == "formal-structure"
    assert stages[3]["status"] == "formal-structure"
    assert stages[7]["status"] == "formal-structure"


def test_endpoint_trace_consistency_remains_non_promoting_and_non_tautological() -> None:
    endpoint = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_section10_endpoint_trace_consistency.json"
    )

    assert endpoint["verified"]["finite_ladder_pointwise_cauchy_consistency"] is True
    assert endpoint["verified"]["majorant_not_refit_from_samples"] is True

    not_verified = endpoint["not_verified"]
    assert all(not_verified.values())
    assert not_verified["actual_section9_residual_jets"] is True
    assert not_verified["uniform_spatial_majorant"] is True
    assert not_verified["locally_uniform_endpoint_limits"] is True
    assert not_verified["all_order_borel_smoothness"] is True
    assert not_verified["smooth_compact_forcing"] is True
    assert not_verified["bounded_energy_and_blowup_closure"] is True
