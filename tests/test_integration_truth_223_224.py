import json
from pathlib import Path

from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "references" / "provenance_manifest_addendum_223_224.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_integration_ledger_maps_prs_223_224_without_promotion() -> None:
    ledger = _load(LEDGER)
    assert ledger["base_main_commit"] == "e109a7aadcc5cbd2cf8c267a672f645ac8e69261"
    assert ledger["integrated_prs"] == [223, 224]
    assert ledger["full_reconstruction"] is False
    assert ledger["paper_exact_velocity_available"] is False

    rows = {row["pr"]: row for row in ledger["layers"]}
    assert set(rows) == {223, 224}
    for row in rows.values():
        assert row["status"] == "formal-structure"
        assert (ROOT / row["artifact"]).is_file()
        assert (ROOT / row["test"]).is_file()
        assert row["capability"].strip()
        assert row["remaining_boundary"].strip()
        if row["provenance"] is not None:
            assert (ROOT / row["provenance"]).is_file()

    assert "fifth-mixed axial U contract" in rows[223]["capability"]
    assert "hierarchy-owned beta fourth-mixed jets" in rows[223]["remaining_boundary"]
    assert "one off-axis point" in rows[224]["capability"]
    assert "non-axisymmetric supplied B jets" in rows[224]["capability"]
    assert rows[224]["provenance"] is None
    assert "truth boundary" in rows[224]["provenance_note"]
    assert "symmetry axis remains uncovered" in rows[224]["remaining_boundary"]
    assert "genuine Navier-Stokes residual/forcing" in rows[224]["remaining_boundary"]


def test_landed_truth_sources_for_223_224_remain_fail_closed() -> None:
    stage2 = _load(ROOT / "references" / "provenance_manifest_addendum_223.json")
    assert stage2["full_reconstruction"] is False
    assert stage2["paper_exact_velocity_available"] is False
    assert stage2["layers"][0]["status"] == "formal-structure"
    assert "AxialFifthMixedJet remains an explicit upstream contract" in stage2["layers"][0]["remaining_boundary"]
    assert "paper-exact velocity remain open" in stage2["layers"][0]["remaining_boundary"]

    regular_flux_source = (
        ROOT
        / "src"
        / "openai_ns_reconstruction"
        / "background_regular_flux_fourth_mixed_jets.py"
    ).read_text(encoding="utf-8")
    assert "class AxialFifthMixedJet:" in regular_flux_source
    assert "The fifth-mixed U jet is an upstream contract" in regular_flux_source
    assert "must not\nbe labelled paper-exact" in regular_flux_source

    divergence_source = (
        ROOT
        / "src"
        / "openai_ns_reconstruction"
        / "section9_section10_divergence_certificate.py"
    ).read_text(encoding="utf-8")
    assert "axis_covered: bool = False" in divergence_source
    assert "actual_section9_sequence_verified: bool = False" in divergence_source
    assert "eq_9_21_infinite_sum_constructed: bool = False" in divergence_source
    assert "section9_field_smooth_extension_through_t1_constructed: bool = False" in divergence_source
    assert "residual_artifact_ready: bool = False" in divergence_source
    assert "endpoint_residual_closure_verified: bool = False" in divergence_source
    assert "finite_energy_verified: bool = False" in divergence_source
    assert "blowup_path_verified: bool = False" in divergence_source
    assert "paper_exact_velocity_available: bool = False" in divergence_source
    assert "axis regularity must come from the actual paper field" in divergence_source
    assert "it is not an analytic-error or convergence bound" in divergence_source

    canonical = _load(ROOT / "references" / "provenance_manifest.json")
    assert canonical["full_reconstruction"] is False

    runtime = construction_status()
    assert runtime["paper_exact_velocity_available"] is False
    assert runtime["status"] == "partial-executable-reconstruction"
    stages = {stage["id"]: stage for stage in runtime["stages"]}
    for stage_id in (1, 2, 3, 6, 7):
        assert stages[stage_id]["status"] == "formal-structure"


def test_cross_stage_forbidden_inferences_and_agent_contract() -> None:
    forbidden = _load(LEDGER)["forbidden_inferences"]
    assert any("fifth-mixed Eq. (5.2)" in x and "Picard convergence" in x for x in forbidden)
    assert any("one-point off-axis" in x and "global divergence-free closure" in x for x in forbidden)
    assert any("symmetry axis" in x and "axis-regularity" in x for x in forbidden)
    assert any("support-exterior zero" in x and "smooth compact forcing" in x for x in forbidden)
    assert any("R-f=0" in x and "independent verification" in x for x in forbidden)
    assert any("Green tests or CI" in x and "paper-exact" in x for x in forbidden)

    agent_contract = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "Never promote status from a green test alone" in agent_contract
    assert "Do not call `f=R(u,p); R-f=0` independent verification" in agent_contract
