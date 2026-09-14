"""Fail-closed bridge to the pinned all-order mixed candidate witness theorem.

This module records a narrower upstream contract than the generic Section 10
endpoint interfaces.  The pinned OpenAI formalization already contains
``MixedCandidateWitness.exists_candidate_witness_of_finite_stages``.  Given the
actual mixed stage families, their genuine ``StageEstimates`` and shrinking
support / one-sided-extension hypotheses, that theorem selects one common
physical schedule and returns the three ``AwayExtensions``, a smooth forcing,
``CandidateProperties``, the downstream consequences, all-order force decay,
and the endpoint force-jet identity.

Python does *not* replay Lean and this file does not manufacture the missing
actual stage data.  It therefore admits only a pinned ``lean-formal-export``
for one exact source/revision and keeps every runtime reconstruction truth flag
false.  Finite endpoint ladders, sampled smoothness, ``f = R`` tautologies, or
formal exports for a different field bundle cannot satisfy this interface.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


PINNED_FORMAL_REPOSITORY = "openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/MixedCandidateWitness.lean"
PINNED_THEOREM = (
    "NavierStokes.MixedCandidateWitness."
    "exists_candidate_witness_of_finite_stages"
)
SECTION10_ENDPOINT_EXACT = Fraction(1, 1)
_UNIVERSAL_ORDER_QUANTIFIER = "forall-natural-derivative-orders"


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class MixedCandidateWitnessFormalExport:
    """Pinned Lean-export witness for the actual mixed Section 9 field bundle.

    Every ``*_certified`` field records an explicit theorem input/output in the
    exported proof object.  None is inferred from numerical evaluation here.
    ``max_derivative_order`` must remain ``None``: the force-decay and endpoint
    jet conclusions of the Lean theorem quantify over every natural order.
    """

    source_id: str
    source_revision: str
    field_bundle_id: str
    stage_estimates_id: str
    potential_stage_family_id: str
    direct_stage_family_id: str
    pressure_stage_family_id: str
    provenance: str
    evidence_kind: str = "lean-formal-export"
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_THEOREM
    endpoint: Fraction = SECTION10_ENDPOINT_EXACT
    derivative_order_quantifier: str = _UNIVERSAL_ORDER_QUANTIFIER
    max_derivative_order: int | None = None
    sampled_or_fitted_evidence: bool = False
    manufactured_residual: bool = False
    same_actual_stage_families_certified: bool = False
    stage_estimates_certified: bool = False
    shrinking_support_all_stages_certified: bool = False
    one_sided_extensions_all_stages_certified: bool = False
    one_common_selected_schedule_certified: bool = False
    selected_schedule_smooth_sums_certified: bool = False
    selected_schedule_vanishing_joint_jets_certified: bool = False
    away_extensions_all_three_sums_certified: bool = False
    smooth_force_constructed_certified: bool = False
    candidate_properties_certified: bool = False
    candidate_consequences_certified: bool = False
    derivative_h3_blowup_certified: bool = False
    force_rapid_decay_all_orders_certified: bool = False
    endpoint_force_jets_all_orders_certified: bool = False

    def __post_init__(self) -> None:
        for name in (
            "source_id",
            "source_revision",
            "field_bundle_id",
            "stage_estimates_id",
            "potential_stage_family_id",
            "direct_stage_family_id",
            "pressure_stage_family_id",
            "provenance",
        ):
            object.__setattr__(self, name, _text(getattr(self, name), name))

        if self.evidence_kind != "lean-formal-export":
            raise ValueError("only lean-formal-export evidence is admissible")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal repository mismatch")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal commit mismatch")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal file mismatch")
        if self.theorem_symbol != PINNED_THEOREM:
            raise ValueError("formal theorem symbol mismatch")
        if type(self.endpoint) is not Fraction or self.endpoint != SECTION10_ENDPOINT_EXACT:
            raise ValueError("endpoint must be the exact rational T=1")
        if self.derivative_order_quantifier != _UNIVERSAL_ORDER_QUANTIFIER:
            raise ValueError("all natural derivative orders must be quantified")
        if self.max_derivative_order is not None:
            raise ValueError("finite derivative-order frontiers are not all-order evidence")
        if self.sampled_or_fitted_evidence is not False:
            raise ValueError("sampled/fitted evidence cannot certify the formal witness")
        if self.manufactured_residual is not False:
            raise ValueError("manufactured residual evidence is forbidden")

        for name in (
            "same_actual_stage_families_certified",
            "stage_estimates_certified",
            "shrinking_support_all_stages_certified",
            "one_sided_extensions_all_stages_certified",
            "one_common_selected_schedule_certified",
            "selected_schedule_smooth_sums_certified",
            "selected_schedule_vanishing_joint_jets_certified",
            "away_extensions_all_three_sums_certified",
            "smooth_force_constructed_certified",
            "candidate_properties_certified",
            "candidate_consequences_certified",
            "derivative_h3_blowup_certified",
            "force_rapid_decay_all_orders_certified",
            "endpoint_force_jets_all_orders_certified",
        ):
            _true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def stage_family_key(self) -> tuple[str, str, str]:
        return (
            self.potential_stage_family_id,
            self.direct_stage_family_id,
            self.pressure_stage_family_id,
        )


@dataclass(frozen=True)
class MixedCandidateWitnessContractCertificate:
    """Accepted formal handoff while runtime truth remains fail-closed."""

    export: MixedCandidateWitnessFormalExport
    expected_source_key: tuple[str, str]
    expected_field_bundle_id: str
    expected_stage_family_key: tuple[str, str, str]
    status: str = "formal-structure"
    actual_section9_sequence_verified: bool = False
    section9_all_order_endpoint_limits_verified: bool = False
    residual_artifact_ready: bool = False
    forcing_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    compact_support_runtime_verified: bool = False
    finite_energy_closure_verified: bool = False
    blow_up_closure_verified: bool = False
    paper_exact_velocity_available: bool = False
    full_reconstruction: bool = False

    @property
    def formal_witness_chain_ready(self) -> bool:
        return (
            self.export.source_key == self.expected_source_key
            and self.export.field_bundle_id == self.expected_field_bundle_id
            and self.export.stage_family_key == self.expected_stage_family_key
            and self.export.theorem_symbol == PINNED_THEOREM
            and self.export.max_derivative_order is None
            and not self.actual_section9_sequence_verified
            and not self.section9_all_order_endpoint_limits_verified
            and not self.residual_artifact_ready
            and not self.forcing_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.compact_support_runtime_verified
            and not self.finite_energy_closure_verified
            and not self.blow_up_closure_verified
            and not self.paper_exact_velocity_available
            and not self.full_reconstruction
        )


def admit_mixed_candidate_witness_export(
    export: MixedCandidateWitnessFormalExport,
    *,
    expected_source_id: str,
    expected_source_revision: str,
    expected_field_bundle_id: str,
    expected_potential_stage_family_id: str,
    expected_direct_stage_family_id: str,
    expected_pressure_stage_family_id: str,
) -> MixedCandidateWitnessContractCertificate:
    """Admit a pinned theorem export only for the already-selected actual data."""

    if not isinstance(export, MixedCandidateWitnessFormalExport):
        raise TypeError("export must be a MixedCandidateWitnessFormalExport")

    expected_source_key = (
        _text(expected_source_id, "expected_source_id"),
        _text(expected_source_revision, "expected_source_revision"),
    )
    expected_bundle = _text(expected_field_bundle_id, "expected_field_bundle_id")
    expected_stage_key = (
        _text(expected_potential_stage_family_id, "expected_potential_stage_family_id"),
        _text(expected_direct_stage_family_id, "expected_direct_stage_family_id"),
        _text(expected_pressure_stage_family_id, "expected_pressure_stage_family_id"),
    )

    if export.source_key != expected_source_key:
        raise ValueError("mixed candidate witness source identity mismatch")
    if export.field_bundle_id != expected_bundle:
        raise ValueError("mixed candidate field bundle identity mismatch")
    if export.stage_family_key != expected_stage_key:
        raise ValueError("mixed candidate stage-family identity mismatch")

    certificate = MixedCandidateWitnessContractCertificate(
        export=export,
        expected_source_key=expected_source_key,
        expected_field_bundle_id=expected_bundle,
        expected_stage_family_key=expected_stage_key,
    )
    if not certificate.formal_witness_chain_ready:
        raise ArithmeticError("mixed candidate witness contract invariant failed")
    return certificate
