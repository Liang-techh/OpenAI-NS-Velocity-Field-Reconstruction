import hashlib
import json
from pathlib import Path

from openai_ns_reconstruction.section9_residual_artifact_fingerprint import (
    Section9ContentAddressedResidualSourceWitness,
    Section9ContentAddressedUniformEnvelopeWitness,
    section9_residual_artifact_sha256,
)
from openai_ns_reconstruction.status import construction_status


ROOT = Path(__file__).resolve().parents[1]


def test_pr165_integration_ledger_preserves_fail_closed_truth() -> None:
    addendum = json.loads(
        (ROOT / "references" / "provenance_manifest_addendum_165_integration.json").read_text(
            encoding="utf-8"
        )
    )
    canonical = json.loads(
        (ROOT / "references" / "provenance_manifest.json").read_text(encoding="utf-8")
    )
    status = construction_status()
    stages = {stage["id"]: stage for stage in status["stages"]}

    assert addendum["base_main_commit"] == "c0e398b4c3c778c5f539b75fd6a064cff9575193"
    assert addendum["integrated_prs"] == [165]
    assert addendum["full_reconstruction"] is False
    assert addendum["paper_exact_velocity_available"] is False
    assert canonical["full_reconstruction"] is False
    assert status["paper_exact_velocity_available"] is False
    assert stages[6]["status"] == "formal-structure"

    layer = addendum["layers"][0]
    assert layer["id"] == "stage-6-residual-improvement-iteration"
    assert layer["status"] == "formal-structure"
    flags = layer["truth_flags"]
    assert flags["preferred_raw_byte_fingerprint_path_available"] is True
    assert flags["caller_supplied_digest_required_on_preferred_path"] is False
    for name in (
        "actual_eq_9_21_residual_verified",
        "source_majorants_derived_from_actual_residual_verified",
        "actual_section9_sequence_verified",
        "proposition_9_9_verified",
        "endpoint_limits_verified",
        "infinite_borel_smoothness_verified",
        "smooth_compact_forcing_verified",
        "bounded_energy_and_blowup_closure_verified",
        "paper_exact_velocity_available",
    ):
        assert flags[name] is False


def test_documented_preferred_path_is_derived_from_exact_artifact_bytes() -> None:
    payload = b"section9-eq921-residual-artifact\nrevision=integration-truth\n"
    mutated = payload[:-2] + b"X\n"

    assert section9_residual_artifact_sha256(payload) == hashlib.sha256(payload).hexdigest()
    assert section9_residual_artifact_sha256(mutated) == hashlib.sha256(mutated).hexdigest()
    assert section9_residual_artifact_sha256(payload) != section9_residual_artifact_sha256(mutated)
    assert callable(Section9ContentAddressedUniformEnvelopeWitness.from_artifact_bytes)
    assert callable(Section9ContentAddressedResidualSourceWitness.from_artifact_bytes)

    provenance = (
        ROOT / "references" / "SECTION9_RESIDUAL_ARTIFACT_FINGERPRINT_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "exact raw serialized residual-artifact bytes" in provenance
    assert "No text decoding, newline normalization, JSON reserialization" in provenance
    assert "audit/content-identity metadata only" in provenance
    assert "does **not** prove that those bytes encode the manuscript's genuine Eq. (9.21) residual" in provenance
    assert "does **not** compute the digest" not in provenance
