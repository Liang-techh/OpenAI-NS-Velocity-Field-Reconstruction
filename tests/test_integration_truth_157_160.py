import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_157_160.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_157_160_and_keeps_truth_gates_closed() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "2f2423fd73727ac6846a5efd32a8d270f86c2f15"
    assert ledger["integrated_prs"] == [157, 158, 159, 160]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {157, 158, 159, 160}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "naturalRemainder(x0)" in rows[157]["remaining_boundary"]
    assert "first genuine Picard iterate" in rows[157]["remaining_boundary"]
    assert "hierarchy-owned partial_eta actualLowerSource/f_n" in rows[158]["remaining_boundary"]
    assert "Proposition 5.3 all-jets residual decay" in rows[158]["remaining_boundary"]
    assert "Proposition 5.5 paper-exact base provider" in rows[159]["remaining_boundary"]
    assert "curl-realized oscillatory correction" in rows[159]["remaining_boundary"]
    assert "audit metadata only" in rows[160]["remaining_boundary"]
    assert "smooth compact forcing" in rows[160]["remaining_boundary"]


def test_source_provenance_canonical_manifest_and_runtime_remain_fail_closed() -> None:
    stage1 = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_natural_operator.json"
    )
    assert stage1["full_reconstruction"] is False
    assert stage1["paper_exact_velocity_available"] is False
    assert stage1["layer"]["status"] == "formal-structure"

    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_158.json")
    assert stage2["full_reconstruction"] is False
    assert stage2["paper_exact_velocity_available"] is False
    assert stage2["layers"][0]["status"] == "formal-structure"

    stage3 = _load(
        ROOT / "references" / "provenance_manifest_addendum_phase_large_band_damping.json"
    )
    assert stage3["full_reconstruction"] is False
    assert stage3["paper_exact_velocity_available"] is False
    assert stage3["layer"]["status"] == "formal-structure"

    stage7_text = (
        ROOT / "references" / "SECTION9_RESIDUAL_ARTIFACT_FINGERPRINT_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "audit metadata only" in stage7_text
    assert "does **not** compute the digest from a manuscript-derived Eq. (9.21) artifact" in stage7_text
    assert "source_majorants_derived_from_actual_residual_verified" in stage7_text
    assert "actual_section9_sequence_verified" in stage7_text
    assert "paper_exact_velocity_available" in stage7_text

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


def test_latest_capabilities_do_not_erase_picard_wave_or_endpoint_boundaries() -> None:
    stage1 = _load(
        ROOT / "references" / "provenance_manifest_addendum_axis_coefficient_natural_operator.json"
    )["layer"]
    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_158.json")["layers"][0]
    stage3 = _load(
        ROOT / "references" / "provenance_manifest_addendum_phase_large_band_damping.json"
    )["layer"]
    stage7_text = (
        ROOT / "references" / "SECTION9_RESIDUAL_ARTIFACT_FINGERPRINT_PROVENANCE.md"
    ).read_text(encoding="utf-8")

    assert "infinite alternating naturalResolvent" in stage1["remaining_boundary"]
    assert "naturalRemainder(x0)" in stage1["remaining_boundary"]
    assert "first genuine Picard iterate" in stage1["remaining_boundary"]

    assert "repaired lower-history hierarchy does not yet own this stronger phi layer" in stage2["remaining_boundary"]
    assert "genuine k=1 Eq. (5.7) Picard application" in stage2["remaining_boundary"]

    assert "No genuine Section 6 carrier/slot family" in stage3["remaining_boundary"]
    assert "actual damping field" in stage3["remaining_boundary"]
    assert "stress-cone wave" in stage3["remaining_boundary"]

    assert "A dishonest caller could still attach the wrong digest" in stage7_text
    assert "actual_section9_sequence_verified" in stage7_text
    assert "endpoint-limit construction/uniqueness" in stage7_text
    assert "infinite Borel right-jet convergence and all-order smoothness" in stage7_text
    assert "smooth compact forcing" in stage7_text
    assert "bounded-energy and blow-up closure" in stage7_text
