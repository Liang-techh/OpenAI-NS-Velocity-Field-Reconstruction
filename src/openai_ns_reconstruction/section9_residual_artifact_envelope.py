"""Load a machine-readable Eq. (9.18) envelope manifest bound to exact residual bytes.

The Section 9 endpoint path already content-addresses a residual artifact and can
analytically convert a uniform Eq. (9.18) envelope into the bounded physical-time
majorant needed near ``t=1``.  Before this adapter, however, callers still had to
construct each ``Section9UniformResidualEnvelopeWitness`` separately from the
artifact bytes.  That leaves a practical integration seam between a materialized
residual artifact and the theorem data claimed for that artifact.

This module closes only that seam.  It consumes two byte strings:

* the exact serialized residual artifact, whose SHA-256 is computed without any
  normalization; and
* a strict JSON theorem-data manifest that names that digest and contains a
  contiguous ``0..N`` family of Eq. (9.18) uniform-envelope data.

The parsed rows are converted to ``Section9UniformResidualEnvelopeWitness``
objects, wrapped with the SHA-256 of the *residual* bytes, and then passed through
the existing analytic endpoint-majorant adapter.  The manifest is not allowed to
replace the residual bytes as the content-addressed object.

This is deliberately not a proof that the bytes are the manuscript's genuine
Eq. (9.21) residual, nor does it derive ``C_{j,m}``, ``K_m``, ``P_{j,m}``, or
flat-remainder estimates from residual algebra.  Those theorem facts remain
external inputs and all paper-exact truth flags stay false.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
import math
from numbers import Integral
from typing import Any

from .section9_endpoint_majorant_adapter import (
    Section9DerivedEndpointMajorantRecord,
    Section9UniformResidualEnvelopeWitness,
    derive_section9_endpoint_majorant_from_uniform_envelope,
)
from .section9_residual_artifact_fingerprint import (
    Section9ContentAddressedUniformEnvelopeWitness,
    section9_residual_artifact_sha256,
)
from .section9_stage_certificate import CertifiedBoundDatum


_SCHEMA = "section9-eq9.21-uniform-envelope-manifest-v1"
_TOP_LEVEL_KEYS = frozenset(
    {
        "schema",
        "residual_artifact_sha256",
        "source_id",
        "source_revision",
        "stage",
        "spatial_window",
        "valid_q_upper",
        "envelopes",
    }
)
_ENVELOPE_KEYS = frozenset(
    {
        "endpoint_derivative_degree",
        "h",
        "derivative_loss",
        "log_power",
        "leading_constant",
        "flat_remainder_uniform_bound",
        "certifications",
    }
)
_BOUND_KEYS = frozenset({"upper_bound", "kind", "provenance"})
_CERTIFICATION_KEYS = frozenset(
    {
        "residual_envelope_uniform_certified",
        "flat_remainder_uniform_certified",
        "official_support_covered_certified",
    }
)


def _strict_object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"duplicate JSON key {key!r} in Section 9 artifact manifest")
        out[key] = value
    return out


def _object(value: object, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"{name} must be a JSON object")
    return value


def _exact_keys(value: dict[str, Any], expected: frozenset[str], name: str) -> None:
    keys = set(value)
    if keys != expected:
        missing = sorted(expected - keys)
        extra = sorted(keys - expected)
        raise ValueError(f"{name} keys mismatch; missing={missing}, extra={extra}")


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _fraction_text(value: object, name: str) -> Fraction:
    if not isinstance(value, str) or not value.strip():
        raise TypeError(f"{name} must be an exact rational string")
    try:
        out = Fraction(value.strip())
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be an exact rational string") from exc
    return out


def _nonnegative_decimal_text(value: object, name: str) -> float:
    if not isinstance(value, str) or not value.strip():
        raise TypeError(f"{name} must be a decimal string")
    try:
        decimal_value = Decimal(value.strip())
    except InvalidOperation as exc:
        raise ValueError(f"{name} must be a finite nonnegative decimal string") from exc
    if not decimal_value.is_finite() or decimal_value < 0:
        raise ValueError(f"{name} must be a finite nonnegative decimal string")
    out = float(decimal_value)
    if not math.isfinite(out):
        raise ValueError(f"{name} is finite in the manifest but not representable in binary64")
    return out


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _bound(value: object, name: str) -> CertifiedBoundDatum:
    row = _object(value, name)
    _exact_keys(row, _BOUND_KEYS, name)
    return CertifiedBoundDatum(
        upper_bound=_nonnegative_decimal_text(row["upper_bound"], f"{name}.upper_bound"),
        kind=_nonempty_text(row["kind"], f"{name}.kind"),
        provenance=_nonempty_text(row["provenance"], f"{name}.provenance"),
    )


@dataclass(frozen=True)
class Section9ResidualArtifactEnvelopeBundle:
    """Parsed theorem data and derived majorants tied to one exact residual artifact."""

    residual_artifact_sha256: str
    theorem_manifest_sha256: str
    source_id: str
    source_revision: str
    stage: int
    spatial_window: int
    envelopes: tuple[Section9UniformResidualEnvelopeWitness, ...]
    content_addressed_envelopes: tuple[
        Section9ContentAddressedUniformEnvelopeWitness, ...
    ]
    derived_majorants: tuple[Section9DerivedEndpointMajorantRecord, ...]
    status: str = "formal-structure"
    manifest_schema_verified: bool = True
    residual_digest_match_verified: bool = True
    contiguous_endpoint_degree_prefix_verified: bool = True
    source_majorants_derived_from_actual_residual_verified: bool = False
    actual_section9_sequence_verified: bool = False
    endpoint_limits_constructed: bool = False
    all_order_borel_smoothness_verified: bool = False
    smooth_compact_forcing_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def max_endpoint_degree(self) -> int:
        return len(self.envelopes) - 1

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def formal_artifact_adapter_ready(self) -> bool:
        degrees = tuple(row.endpoint_derivative_degree for row in self.envelopes)
        return (
            bool(self.envelopes)
            and degrees == tuple(range(len(self.envelopes)))
            and len(self.envelopes)
            == len(self.content_addressed_envelopes)
            == len(self.derived_majorants)
            and all(row.source_key == self.source_key for row in self.envelopes)
            and all(row.stage == self.stage for row in self.envelopes)
            and all(row.spatial_window == self.spatial_window for row in self.envelopes)
            and all(
                row.residual_artifact_sha256 == self.residual_artifact_sha256
                for row in self.content_addressed_envelopes
            )
            and all(row.formal_adapter_ready for row in self.derived_majorants)
            and self.manifest_schema_verified
            and self.residual_digest_match_verified
            and self.contiguous_endpoint_degree_prefix_verified
            and not self.source_majorants_derived_from_actual_residual_verified
            and not self.actual_section9_sequence_verified
            and not self.endpoint_limits_constructed
            and not self.all_order_borel_smoothness_verified
            and not self.smooth_compact_forcing_verified
            and not self.paper_exact_velocity_available
        )


def derive_section9_endpoint_majorants_from_residual_artifact(
    residual_artifact: bytes,
    theorem_manifest: bytes,
) -> Section9ResidualArtifactEnvelopeBundle:
    """Parse one strict theorem manifest and derive its finite endpoint ladder.

    ``residual_artifact`` is the content-addressed object.  ``theorem_manifest``
    is a sidecar whose declared digest must match those exact bytes.  The JSON
    parser rejects duplicate keys and every object uses an exact key set so an
    unrecognized field cannot silently alter the artifact interpretation.

    The manifest must carry a contiguous endpoint-degree prefix ``0..N``.  Each
    row is fed through the existing analytic Eq. (9.18) -> endpoint-majorant
    adapter, so the output coefficient is machine-derived from the manifest's
    theorem data rather than independently supplied.
    """

    residual_digest = section9_residual_artifact_sha256(residual_artifact)
    if not isinstance(theorem_manifest, bytes):
        raise TypeError("theorem_manifest must be exact bytes")
    if not theorem_manifest:
        raise ValueError("theorem_manifest must be nonempty")
    theorem_manifest_sha256 = hashlib.sha256(theorem_manifest).hexdigest()

    try:
        text = theorem_manifest.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("theorem_manifest must be UTF-8 JSON") from exc
    try:
        payload = json.loads(text, object_pairs_hook=_strict_object_pairs)
    except json.JSONDecodeError as exc:
        raise ValueError("theorem_manifest must be valid JSON") from exc

    top = _object(payload, "theorem manifest")
    _exact_keys(top, _TOP_LEVEL_KEYS, "theorem manifest")
    if top["schema"] != _SCHEMA:
        raise ValueError(f"theorem manifest schema must be {_SCHEMA!r}")
    declared_digest = _nonempty_text(
        top["residual_artifact_sha256"], "residual_artifact_sha256"
    )
    if declared_digest != residual_digest:
        raise ValueError(
            "theorem manifest residual_artifact_sha256 does not match exact residual bytes"
        )

    source_id = _nonempty_text(top["source_id"], "source_id")
    source_revision = _nonempty_text(top["source_revision"], "source_revision")
    stage = _natural(top["stage"], "stage")
    spatial_window = _natural(top["spatial_window"], "spatial_window")
    valid_q_upper = _fraction_text(top["valid_q_upper"], "valid_q_upper")

    raw_envelopes = top["envelopes"]
    if not isinstance(raw_envelopes, list) or not raw_envelopes:
        raise ValueError("envelopes must be a nonempty JSON array")

    by_degree: dict[int, Section9UniformResidualEnvelopeWitness] = {}
    for index, raw in enumerate(raw_envelopes):
        row = _object(raw, f"envelopes[{index}]")
        _exact_keys(row, _ENVELOPE_KEYS, f"envelopes[{index}]")
        degree = _natural(
            row["endpoint_derivative_degree"],
            f"envelopes[{index}].endpoint_derivative_degree",
        )
        if degree in by_degree:
            raise ValueError(f"duplicate endpoint_derivative_degree {degree}")

        certifications = _object(
            row["certifications"], f"envelopes[{index}].certifications"
        )
        _exact_keys(
            certifications,
            _CERTIFICATION_KEYS,
            f"envelopes[{index}].certifications",
        )

        witness = Section9UniformResidualEnvelopeWitness(
            endpoint_derivative_degree=degree,
            stage=stage,
            spatial_window=spatial_window,
            h=_fraction_text(row["h"], f"envelopes[{index}].h"),
            derivative_loss=_fraction_text(
                row["derivative_loss"], f"envelopes[{index}].derivative_loss"
            ),
            log_power=_natural(row["log_power"], f"envelopes[{index}].log_power"),
            leading_constant=_bound(
                row["leading_constant"], f"envelopes[{index}].leading_constant"
            ),
            flat_remainder_uniform_bound=_bound(
                row["flat_remainder_uniform_bound"],
                f"envelopes[{index}].flat_remainder_uniform_bound",
            ),
            valid_q_upper=valid_q_upper,
            source_id=source_id,
            source_revision=source_revision,
            residual_envelope_uniform_certified=_strict_true(
                certifications["residual_envelope_uniform_certified"],
                f"envelopes[{index}].residual_envelope_uniform_certified",
            ),
            flat_remainder_uniform_certified=_strict_true(
                certifications["flat_remainder_uniform_certified"],
                f"envelopes[{index}].flat_remainder_uniform_certified",
            ),
            official_support_covered_certified=_strict_true(
                certifications["official_support_covered_certified"],
                f"envelopes[{index}].official_support_covered_certified",
            ),
        )
        by_degree[degree] = witness

    required = set(range(len(by_degree)))
    supplied = set(by_degree)
    if supplied != required:
        raise ValueError(
            "artifact envelope ladder must cover a contiguous endpoint-degree prefix 0..N; "
            f"missing={sorted(required - supplied)}, extra={sorted(supplied - required)}"
        )
    envelopes = tuple(by_degree[degree] for degree in range(len(by_degree)))
    content_addressed = tuple(
        Section9ContentAddressedUniformEnvelopeWitness.from_artifact_bytes(
            envelope=row,
            residual_artifact=residual_artifact,
        )
        for row in envelopes
    )
    derived = tuple(
        derive_section9_endpoint_majorant_from_uniform_envelope(row)
        for row in envelopes
    )

    result = Section9ResidualArtifactEnvelopeBundle(
        residual_artifact_sha256=residual_digest,
        theorem_manifest_sha256=theorem_manifest_sha256,
        source_id=source_id,
        source_revision=source_revision,
        stage=stage,
        spatial_window=spatial_window,
        envelopes=envelopes,
        content_addressed_envelopes=content_addressed,
        derived_majorants=derived,
    )
    if not result.formal_artifact_adapter_ready:
        raise ArithmeticError("Section 9 residual-artifact envelope adapter invariant failed")
    return result
