import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_188_191.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_188_191_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "fa100ccc5ebbea31996dd860ef4500ca84e775b5"
    assert ledger["integrated_prs"] == [188, 189, 190, 191]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {188, 189, 190, 191}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert (ROOT / row["provenance"]).is_file()
        assert row["equation_or_source"].strip()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()

    assert "mixed-scale decomposition" in rows[188]["capability"]
    assert "naturalRemainder(x0)" in rows[188]["remaining_boundary"]
    assert "partial_eta A0" in rows[189]["capability"]
    assert "partial_eta^2 W_n^(0)" in rows[189]["remaining_boundary"]
    assert "exact same active scope witness" in rows[190]["capability"]
    assert "CellIndex" in rows[190]["remaining_boundary"]
    assert "multivariate Taylor composition" in rows[191]["capability"]
    assert "provider inputs" in rows[191]["remaining_boundary"]


def test_source_provenance_and_runtime_remain_fail_closed() -> None:
    axial = _load(
        ROOT
        / "references"
        / "provenance_manifest_addendum_axis_coefficient_wide_axial_remainder.json"
    )
    assert axial["full_reconstruction"] is False
    assert axial["paper_exact_velocity_available"] is False
    assert axial["layer"]["status"] == "formal-structure"
    assert axial["layer"]["actual_data"]["axial_formula"] == (
        "inverseL * (lin2 - (1/Lambda) * slow2 + pressure)"
    )
    assert "angular naturalRemainder branch" in axial["layer"]["remaining_boundary"]
    assert "Picard x1" in axial["layer"]["remaining_boundary"]

    matrix_jets = _load(ROOT / "references" / "provenance_manifest_addendum_189.json")
    assert matrix_jets["full_reconstruction"] is False
    assert matrix_jets["paper_exact_velocity_available"] is False
    assert matrix_jets["layers"][0]["status"] == "formal-structure"
    assert "partial_eta A0" in matrix_jets["layers"][0]["landed_capability"]
    assert "partial_eta^2 W_n^(0)" in matrix_jets["layers"][0]["remaining_boundary"]

    physical_estimate = _load(
        ROOT / "references" / "PHASE_LARGE_BAND_PHYSICAL_ESTIMATE_COVERAGE_PROVENANCE.json"
    )
    assert physical_estimate["full_reconstruction"] is False
    assert physical_estimate["paper_exact_velocity_available"] is False
    assert physical_estimate["layer"]["status"] == "formal-structure"
    truth = physical_estimate["layer"]["truth_boundary"]
    assert truth["physical_estimate_scope_machine_checked"] is True
    assert truth["actual_cell_index_enumeration_machine_verified"] is False
    assert truth["actual_physical_base_values_materialized"] is False
    assert truth["theorem_applications_machine_replayed"] is False
    assert truth["uniform_eq_7_9_to_7_11_verified"] is False
    assert truth["paper_exact_velocity_available"] is False

    cutoff_provenance = (
        ROOT / "references" / "SECTION9_CUTOFF_COMPOSITION_JET_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "similarity_coordinate_jet_machine_derived_from_eq_4_1 = false" in cutoff_provenance
    assert "paper_fixed_cutoff_derivatives_machine_verified = false" in cutoff_provenance
    assert "paper_exact_velocity_available = false" in cutoff_provenance
    assert "q-jet" in cutoff_provenance and "provider inputs" in cutoff_provenance

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_are_explicit() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]

    assert any(
        "axial" in item and "naturalRemainder(x0)" in item and "fixed point" in item
        for item in forbidden
    )
    assert any(
        "A0/A1" in item and "Picard convergence" in item and "Proposition 5.3" in item
        for item in forbidden
    )
    assert any(
        "CellIndex" in item and "FinalSlowBase.velocity" in item and "Eqs. (7.9)-(7.11)" in item
        for item in forbidden
    )
    assert any(
        "provider-supplied" in item
        and "infinite Eq. (9.21)" in item
        and "residual artifact" in item
        for item in forbidden
    )
    assert any("R-f=0" in item and "independent verification" in item for item in forbidden)
    assert any("Green tests or CI" in item and "paper-exact" in item for item in forbidden)
