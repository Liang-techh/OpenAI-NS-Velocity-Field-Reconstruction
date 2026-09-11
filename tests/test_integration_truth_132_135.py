import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_132_135.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_landed_artifacts_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "a34622932b3bd893acd78c8b503605f9448f828b"
    assert ledger["integrated_prs"] == [132, 133, 134, 135]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {132, 133, 134, 135}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["remaining_boundary"].strip()

    assert "complete coefficientOperators record" in rows[132]["remaining_boundary"]
    assert "naturalRemainder(x0)" in rows[132]["remaining_boundary"]
    assert "Lemma 5.2 repair must be lifted" in rows[133]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in rows[133]["remaining_boundary"]
    assert "No paper-exact Proposition 5.5 background" in rows[134]["remaining_boundary"]
    assert "Actual Section 9 residual derivative limits" in rows[135]["remaining_boundary"]
    assert "infinite smooth extension" in rows[135]["remaining_boundary"]


def test_source_addenda_canonical_manifest_and_runtime_remain_fail_closed() -> None:
    addenda = (
        "provenance_manifest_addendum_axis_coefficient_inverse_param_product.json",
        "provenance_manifest_addendum_133.json",
        "provenance_manifest_addendum_phase_large_band_local_base.json",
        "provenance_manifest_addendum_section10_endpoint_borel_admission.json",
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


def test_section10_endpoint_candidate_admission_remains_strictly_non_promoting() -> None:
    endpoint = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_section10_endpoint_borel_admission.json"
    )
    layer = endpoint["layer"]

    verified = layer["verified"]
    assert all(verified.values())
    assert verified["finite_endpoint_trace_enclosure_membership_gate"] is True
    assert verified["candidate_frozen_between_left_and_right_checks"] is True
    assert verified["rejection_occurs_before_borel_scale_query"] is True
    assert verified["finite_prefix_right_jet_admission_after_trace_consistency"] is True

    not_verified = layer["not_verified"]
    assert all(not_verified.values())
    assert not_verified["source_majorants_derived_from_actual_section9_residual"] is True
    assert not_verified["actual_section9_residual_endpoint_limits"] is True
    assert not_verified["endpoint_limit_uniqueness"] is True
    assert not_verified["analytic_template_bounds"] is True
    assert not_verified["infinite_borel_right_jets"] is True
    assert not_verified["all_order_borel_smoothness"] is True
    assert not_verified["smooth_compact_forcing"] is True
    assert not_verified["bounded_energy_and_blowup_closure"] is True
