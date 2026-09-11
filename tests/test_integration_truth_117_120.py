import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_117_120.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_landed_artifacts_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "0ebce7b28faac20d14399915c282a6c5f0bb7dcd"
    assert ledger["integrated_prs"] == [117, 118, 119, 120]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {117, 118, 119, 120}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["remaining_boundary"].strip()

    assert "naturalRemainder" in rows[117]["remaining_boundary"]
    assert "Only the k=0 Picard term" in rows[118]["remaining_boundary"]
    assert "no infinite Eq. (9.21) field" in rows[119]["remaining_boundary"]
    assert "does not cover t=1" in rows[120]["remaining_boundary"]


def test_source_addenda_and_runtime_remain_fail_closed() -> None:
    addenda = (
        "provenance_manifest_addendum_axis_coefficient_parameter_primitive.json",
        "provenance_manifest_addendum_118.json",
        "provenance_manifest_addendum_section9_finite_prefix.json",
        "provenance_manifest_addendum_section9_physical_prefix_completion.json",
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


def test_section9_prefix_layers_preserve_non_tautological_boundaries() -> None:
    finite = _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_finite_prefix.json"
    )
    physical = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_section9_physical_prefix_completion.json"
    )

    finite_boundary = finite["layer"]["remaining_boundary"]
    assert "Section 7/8 correction values" in finite_boundary
    assert "external cutoff callable" in finite_boundary
    assert "finite prefix" in finite_boundary
    assert "infinite Eq. (9.21) locally finite sum" in finite_boundary

    not_verified = physical["not_verified"]
    assert all(not_verified.values())
    assert not_verified["actual_correction_field_values"] is True
    assert not_verified["paper_fixed_cutoff_pointwise_evaluator"] is True
    assert not_verified["eq_9_21_infinite_field_constructed"] is True
    assert not_verified["endpoint_t_equals_1_covered"] is True
    assert not_verified["endpoint_uniform_residual_majorants"] is True
    assert not_verified["smooth_compact_forcing"] is True
    assert not_verified["bounded_energy_and_blowup_closure"] is True
