import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_127_130.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_landed_artifacts_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "b75571dd4995cdfd2ffa840c1416e69ab8ea4594"
    assert ledger["integrated_prs"] == [127, 128, 129, 130]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {127, 128, 129, 130}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["remaining_boundary"].strip()

    assert "complete coefficientOperators" in rows[127]["remaining_boundary"]
    assert "naturalRemainder(x0)" in rows[127]["remaining_boundary"]
    assert "Eq. (5.2) must be lifted" in rows[128]["remaining_boundary"]
    assert "No paper-exact Proposition 5.5 background" in rows[129]["remaining_boundary"]
    assert "infinite rightExtension_right_jets theorem" in rows[130]["remaining_boundary"]


def test_source_addenda_canonical_manifest_and_runtime_remain_fail_closed() -> None:
    addenda = (
        "provenance_manifest_addendum_axis_coefficient_regular_inverse.json",
        "provenance_manifest_addendum_128.json",
        "provenance_manifest_addendum_phase_large_band_frame.json",
        "provenance_manifest_addendum_section10_borel_prefix_right_jets.json",
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


def test_section10_finite_prefix_right_jets_remain_strictly_non_promoting() -> None:
    endpoint = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_section10_borel_prefix_right_jets.json"
    )
    layer = endpoint["layer"]

    verified = layer["verified"]
    assert verified["dense_full_spacetime_jet_shape_gate"] is True
    assert verified["finite_prefix_common_unit_plateau"] is True
    assert verified["finite_prefix_endpoint_normal_jet_algebra"] is True
    assert verified["independent_future_branch_interpolation_regression"] is True

    not_verified = layer["not_verified"]
    assert all(not_verified.values())
    assert not_verified["actual_section9_residual_endpoint_limits"] is True
    assert not_verified["analytic_template_bounds"] is True
    assert not_verified["infinite_borel_right_jets"] is True
    assert not_verified["all_order_borel_smoothness"] is True
    assert not_verified["smooth_compact_forcing"] is True
    assert not_verified["bounded_energy_and_blowup_closure"] is True
