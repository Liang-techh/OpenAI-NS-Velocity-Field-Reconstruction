from fractions import Fraction
import hashlib
import json

import numpy as np
import pytest

from openai_ns_reconstruction.section9_eq921_prefix_source_artifact import (
    materialize_eq_9_21_finite_prefix_source_artifact,
)
from openai_ns_reconstruction.section9_finite_prefix import (
    Section9FinitePrefixEvaluationCertificate,
    Section9FinitePrefixStageContribution,
)


def _readonly_vec(values):
    out = np.asarray(values, dtype=float)
    out.setflags(write=False)
    return out


def _prefix(*, stage2_b_term=8.0, actual_values_verified=False):
    contributions = (
        Section9FinitePrefixStageContribution(
            stage=1,
            scale=8,
            cutoff_argument=Fraction(1, 5),
            cutoff_weight=1.0,
            A_term=_readonly_vec([1.0, 2.0, -1.0]),
            B_term=4.0,
            p_term=-0.5,
            value_provenance="test stage 1; not paper correction data",
            support_zero_short_circuit=False,
        ),
        Section9FinitePrefixStageContribution(
            stage=2,
            scale=16,
            cutoff_argument=Fraction(2, 5),
            cutoff_weight=0.5,
            A_term=_readonly_vec([2.0, -3.0, 0.25]),
            B_term=stage2_b_term,
            p_term=1.5,
            value_provenance="test stage 2; not paper correction data",
            support_zero_short_circuit=False,
        ),
    )
    return Section9FinitePrefixEvaluationCertificate(
        q=Fraction(1, 40),
        q_big=Fraction(1, 4),
        prefix_order=2,
        stages=(1, 2),
        contributions=contributions,
        A_correction_prefix=_readonly_vec([3.0, -1.0, -0.75]),
        B_correction_prefix=4.0 + stage2_b_term,
        p_correction_prefix=1.0,
        cutoff_evaluator_provenance="test-only cutoff; not manuscript pointwise data",
        actual_correction_field_values_verified=actual_values_verified,
    )


def test_materializes_deterministic_content_addressed_source_bytes():
    first = materialize_eq_9_21_finite_prefix_source_artifact(_prefix())
    second = materialize_eq_9_21_finite_prefix_source_artifact(_prefix())

    assert first.artifact_bytes == second.artifact_bytes
    assert first.sha256 == hashlib.sha256(first.artifact_bytes).hexdigest()
    assert first.sha256 == second.sha256
    assert first.formal_source_artifact_ready

    payload = json.loads(first.artifact_bytes)
    assert payload["schema"] == "section9-eq9.21-finite-prefix-source-artifact-v1"
    assert payload["artifact_kind"] == "eq9.21-finite-prefix-source-not-residual"
    assert payload["q"] == {"numerator": 1, "denominator": 40}
    assert payload["q_big"] == {"numerator": 1, "denominator": 4}
    assert payload["contributions"][0]["cutoff_argument"] == {
        "numerator": 1,
        "denominator": 5,
    }
    assert payload["contributions"][0]["cutoff_weight_binary64"] == float(1.0).hex()
    assert payload["truth_boundary"]["residual_artifact_ready"] is False
    assert payload["truth_boundary"]["paper_exact_velocity_available"] is False

    assert not first.residual_artifact_ready
    assert not first.source_majorants_derived_from_actual_residual_verified
    assert not first.paper_exact_velocity_available


def test_payload_change_changes_source_artifact_identity():
    baseline = materialize_eq_9_21_finite_prefix_source_artifact(_prefix())
    changed = materialize_eq_9_21_finite_prefix_source_artifact(
        _prefix(stage2_b_term=8.25)
    )

    assert baseline.artifact_bytes != changed.artifact_bytes
    assert baseline.sha256 != changed.sha256


def test_binary64_payload_is_canonical_not_json_decimal_dependent():
    artifact = materialize_eq_9_21_finite_prefix_source_artifact(_prefix())
    payload = json.loads(artifact.artifact_bytes)

    assert payload["A_correction_prefix_binary64"] == [
        float(3.0).hex(),
        float(-1.0).hex(),
        float(-0.75).hex(),
    ]
    assert payload["contributions"][1]["B_term_binary64"] == float(8.0).hex()


def test_rejects_prefix_that_claims_actual_correction_values_verified():
    prefix = _prefix(actual_values_verified=True)
    assert not prefix.formal_prefix_ready
    with pytest.raises(ValueError, match="not formally ready"):
        materialize_eq_9_21_finite_prefix_source_artifact(prefix)


def test_rejects_non_prefix_input():
    with pytest.raises(
        TypeError, match="Section9FinitePrefixEvaluationCertificate"
    ):
        materialize_eq_9_21_finite_prefix_source_artifact(object())
