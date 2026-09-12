"""Canonical artifact materialization for an admitted finite Eq. (9.21) prefix.

This is a provider-side step toward a genuine Section 9 residual artifact.  It
takes the already-landed :class:`Section9FinitePrefixEvaluationCertificate` and
serializes its exact finite-prefix source data into deterministic UTF-8 JSON
bytes.  Exact rationals are written as numerator/denominator pairs and binary64
values are written with ``float.hex()`` so content identity does not depend on
JSON floating-point formatting.

The resulting bytes are deliberately labelled *source*, not *residual*.  A
finite pointwise prefix contains no spacetime derivative data and therefore
cannot by itself define the Navier--Stokes residual or the Eq. (9.18) theorem
constants.  In particular, these bytes must not be passed off as the genuine
Eq. (9.21) residual artifact accepted by
``section9_residual_artifact_envelope``.  All paper-exact gates remain false.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from typing import Any

import numpy as np

from .section9_finite_prefix import Section9FinitePrefixEvaluationCertificate


_SCHEMA = "section9-eq9.21-finite-prefix-source-artifact-v1"
_ARTIFACT_KIND = "eq9.21-finite-prefix-source-not-residual"


def _fraction_payload(value: Fraction) -> dict[str, int]:
    if not isinstance(value, Fraction):
        raise TypeError("finite-prefix rational payload must be Fraction")
    return {"numerator": value.numerator, "denominator": value.denominator}


def _float_hex(value: object, name: str) -> str:
    out = float(value)
    if not np.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out.hex()


def _vector_hex(value: object, name: str) -> list[str]:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    return [float(component).hex() for component in out]


def _canonical_json_bytes(payload: dict[str, Any]) -> bytes:
    text = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return text.encode("utf-8")


@dataclass(frozen=True)
class Section9Eq921FinitePrefixSourceArtifact:
    """Content-addressed finite-prefix *source* artifact.

    ``residual_artifact_ready`` is permanently false because the input
    certificate carries only pointwise finite-prefix values.  The artifact is
    useful as a stable upstream identity for a future analytic residual
    provider, but it cannot substitute for that provider.
    """

    artifact_bytes: bytes
    sha256: str
    prefix_order: int
    q: Fraction
    q_big: Fraction
    artifact_kind: str = _ARTIFACT_KIND
    status: str = "formal-structure"
    finite_prefix_source_materialized: bool = True
    residual_artifact_ready: bool = False
    actual_correction_field_values_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    source_majorants_derived_from_actual_residual_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_source_artifact_ready(self) -> bool:
        return (
            bool(self.artifact_bytes)
            and len(self.sha256) == 64
            and self.artifact_kind == _ARTIFACT_KIND
            and self.status == "formal-structure"
            and self.finite_prefix_source_materialized
            and not self.residual_artifact_ready
            and not self.actual_correction_field_values_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.source_majorants_derived_from_actual_residual_verified
            and not self.paper_exact_velocity_available
        )


def materialize_eq_9_21_finite_prefix_source_artifact(
    prefix: Section9FinitePrefixEvaluationCertificate,
) -> Section9Eq921FinitePrefixSourceArtifact:
    """Serialize one admitted finite Eq. (9.21) prefix deterministically.

    Only a fail-closed ``formal_prefix_ready`` certificate is accepted.  The
    serializer does not call a residual evaluator, does not differentiate the
    fields, and does not consume an Eq. (9.18) sidecar.
    """

    if not isinstance(prefix, Section9FinitePrefixEvaluationCertificate):
        raise TypeError(
            "prefix must be a Section9FinitePrefixEvaluationCertificate"
        )
    if not prefix.formal_prefix_ready:
        raise ValueError("finite Eq. (9.21) prefix is not formally ready")

    rows: list[dict[str, Any]] = []
    for contribution in prefix.contributions:
        rows.append(
            {
                "stage": contribution.stage,
                "scale": contribution.scale,
                "cutoff_argument": _fraction_payload(contribution.cutoff_argument),
                "cutoff_weight_binary64": _float_hex(
                    contribution.cutoff_weight, "cutoff_weight"
                ),
                "A_term_binary64": _vector_hex(contribution.A_term, "A_term"),
                "B_term_binary64": _float_hex(contribution.B_term, "B_term"),
                "p_term_binary64": _float_hex(contribution.p_term, "p_term"),
                "value_provenance": contribution.value_provenance,
                "support_zero_short_circuit": bool(
                    contribution.support_zero_short_circuit
                ),
            }
        )

    payload = {
        "schema": _SCHEMA,
        "artifact_kind": _ARTIFACT_KIND,
        "status": "formal-structure",
        "q": _fraction_payload(prefix.q),
        "q_big": _fraction_payload(prefix.q_big),
        "prefix_order": prefix.prefix_order,
        "stages": list(prefix.stages),
        "cutoff_evaluator_provenance": prefix.cutoff_evaluator_provenance,
        "contributions": rows,
        "A_correction_prefix_binary64": _vector_hex(
            prefix.A_correction_prefix, "A_correction_prefix"
        ),
        "B_correction_prefix_binary64": _float_hex(
            prefix.B_correction_prefix, "B_correction_prefix"
        ),
        "p_correction_prefix_binary64": _float_hex(
            prefix.p_correction_prefix, "p_correction_prefix"
        ),
        "truth_boundary": {
            "actual_correction_field_values_verified": False,
            "paper_fixed_cutoff_pointwise_evaluator_verified": False,
            "infinite_schedule_constructed_here": False,
            "eq_9_21_infinite_sum_constructed": False,
            "residual_artifact_ready": False,
            "source_majorants_derived_from_actual_residual_verified": False,
            "paper_exact_velocity_available": False,
        },
    }
    artifact_bytes = _canonical_json_bytes(payload)
    result = Section9Eq921FinitePrefixSourceArtifact(
        artifact_bytes=artifact_bytes,
        sha256=hashlib.sha256(artifact_bytes).hexdigest(),
        prefix_order=prefix.prefix_order,
        q=prefix.q,
        q_big=prefix.q_big,
    )
    if not result.formal_source_artifact_ready:
        raise ArithmeticError("finite-prefix source-artifact invariant failed")
    return result
