import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_170_171.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_170_171_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "7b41bc43d09e5dcc8fdbc691117a17dac5d95d4c"
    assert ledger["integrated_prs"] == [170, 171]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {170, 171}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "CellIndex" in rows[170]["remaining_boundary"]
    assert "Proposition 5.5" in rows[170]["remaining_boundary"]
    assert "curl-realized oscillatory waves" in rows[170]["remaining_boundary"]
    assert "genuine Eq. (9.21) residual" in rows[171]["remaining_boundary"]
    assert "all-order endpoint limits" in rows[171]["remaining_boundary"]
    assert "smooth compact forcing" in rows[171]["remaining_boundary"]


def test_source_provenance_canonical_manifest_and_runtime_stay_fail_closed() -> None:
    primary_geometry = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_PRIMARY_GEOMETRY_SCOPE_PROVENANCE.json"
    )
    assert primary_geometry["status"] == "formal-structure"
    flags = primary_geometry["truth_flags"]
    assert flags["actual_cell_index_enumeration_machine_verified"] is False
    assert flags["actual_active_family_scope_machine_verified"] is False
    assert flags["actual_base_fields_verified"] is False
    assert flags["uniform_eq_7_9_to_7_11_verified"] is False
    assert flags["paper_exact_velocity_available"] is False

    residual_envelope = _load(
        ROOT / "references" / "provenance_manifest_addendum_section9_artifact_envelope_adapter.json"
    )
    assert residual_envelope["status"] == "formal-structure"
    assert residual_envelope["source_majorants_derived_from_actual_residual_verified"] is False
    assert residual_envelope["actual_section9_sequence_verified"] is False
    assert residual_envelope["endpoint_limits_constructed"] is False
    assert residual_envelope["all_order_borel_smoothness_verified"] is False
    assert residual_envelope["smooth_compact_forcing_verified"] is False
    assert residual_envelope["full_reconstruction"] is False
    assert residual_envelope["paper_exact_velocity_available"] is False

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    assert stages[3]["status"] == "formal-structure"
    assert stages[6]["status"] == "formal-structure"
    assert stages[7]["status"] == "formal-structure"


def test_new_bindings_do_not_close_scope_endpoint_or_verification_gates() -> None:
    ledger = _load(LEDGER)
    rows = {row["pr"]: row for row in ledger["layers"]}

    assert "both signs are reconstructed mechanically" in rows[170]["capability"]
    assert "does not enumerate" in rows[170]["remaining_boundary"]
    assert "actual manuscript active family" in rows[170]["remaining_boundary"]

    assert "Exact residual artifact bytes" in rows[171]["capability"]
    assert "sidecar still supplies C_{j,m}, K_m, P_{j,m}" in rows[171]["remaining_boundary"]
    assert "source majorants derived from the actual residual" in rows[171]["remaining_boundary"]

    forbidden = ledger["forbidden_inferences"]
    assert any("CellIndex" in item and "machine verified" in item for item in forbidden)
    assert any("SHA-256" in item and "Eq. (9.21) residual" in item for item in forbidden)
    assert any("finite contiguous endpoint-majorant ladder" in item for item in forbidden)
    assert any("setting f equal to the evaluated residual" in item for item in forbidden)
