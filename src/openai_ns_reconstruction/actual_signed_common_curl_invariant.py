"""Bind the actual common-curl theorem to jets derived from the cycle invariant.

``ActualSignedCommonDynamics.common_curl_and_divergence`` consumes a
``UniformLocalJets`` hypothesis for an otherwise generic request.  The pinned
formal source also proves ``ActualSignedStageControls.fullRequest_jets_from_invariant``:
the actual cycle analytic invariant supplies exactly those jets for
``LocalSignedRequest.fullRequest``.

This module is a fail-closed *composition* gate.  It does not replay Lean and
does not materialize any wave.  It records that the jets theorem is instantiated
with ``K = phaseCell`` and with the same request, exponent, and ``B/N0`` used by
an already admitted common-curl application.  Sampled, fitted, or numerical
metadata cannot be promoted through this interface.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .actual_signed_common_curl import (
    ActualSignedCommonCurlAdmission,
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
)

PINNED_JETS_FILE = "NavierStokes/ActualSignedStageControls.lean"
PINNED_JETS_THEOREM = "ActualSignedStageControls.fullRequest_jets_from_invariant"
PINNED_FULL_REQUEST = "LocalSignedRequest.fullRequest"
PINNED_PHASE_CELL = "ActualSignedStageControls.phaseCell"
PINNED_FULL_STRIP = "ActualSignedStageControls.fullStrip"
PINNED_ACTUAL_STRIP = "ActualPrimaryBounds.strip"
PINNED_CYCLE_INVARIANT = "CorrectionStep.CycleAnalyticInvariant"
PINNED_UNIFORM_LOCAL_JETS = "PeriodizedWaveBounds.UniformLocalJets"
PINNED_PRODUCT_STRIP = "HarmonicWaveInteraction.productStrip"

REQUIRED_JETS_SYMBOLS = (
    PINNED_FULL_REQUEST,
    PINNED_PHASE_CELL,
    PINNED_FULL_STRIP,
    PINNED_ACTUAL_STRIP,
    PINNED_CYCLE_INVARIANT,
    PINNED_UNIFORM_LOCAL_JETS,
    PINNED_PRODUCT_STRIP,
)

_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _pinned(value: object, expected: str, name: str) -> str:
    text = _nonempty_text(value, name)
    if text != expected:
        raise ValueError(f"{name} must match the pinned formal source exactly")
    return text


@dataclass(frozen=True)
class ActualSignedInvariantJetsWitness:
    """One formal application of ``fullRequest_jets_from_invariant``.

    Lean objects are deliberately represented by opaque stable strings.  Python
    never tries to reconstruct a ``CycleState`` or ``SignedLabel``.  The
    boolean identity facts are expected to come from theorem/export metadata,
    not from numerical comparison.
    """

    B: int
    N0: int
    request_repr: str
    beta_repr: str
    cycle_state_repr: str
    geometry_repr: str
    context_repr: str
    primary_repr: str
    pressure_repr: str
    label_carrier_repr: str
    application_id: str
    evidence_kind: str
    provenance: str
    cycle_analytic_invariant_certified: bool
    full_request_identity_certified: bool
    geometry_strip_is_actual_strip_certified: bool
    phase_cell_instantiation_certified: bool
    beta_is_invariant_sigma_certified: bool
    theorem_application_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_JETS_FILE
    theorem: str = PINNED_JETS_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_JETS_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        for name in (
            "request_repr",
            "beta_repr",
            "cycle_state_repr",
            "geometry_repr",
            "context_repr",
            "primary_repr",
            "pressure_repr",
            "label_carrier_repr",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )

        for name in (
            "cycle_analytic_invariant_certified",
            "full_request_identity_certified",
            "geometry_strip_is_actual_strip_certified",
            "phase_cell_instantiation_certified",
            "beta_is_invariant_sigma_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)

        for name, expected in (
            ("lean_repository", PINNED_LEAN_REPOSITORY),
            ("lean_commit", PINNED_LEAN_COMMIT),
            ("lean_file", PINNED_JETS_FILE),
            ("theorem", PINNED_JETS_THEOREM),
        ):
            object.__setattr__(self, name, _pinned(getattr(self, name), expected, name))

        if not isinstance(self.dependency_symbols, tuple):
            raise ValueError("dependency_symbols must be the pinned theorem-symbol tuple")
        symbols = tuple(
            _nonempty_text(symbol, "dependency symbol") for symbol in self.dependency_symbols
        )
        if symbols != REQUIRED_JETS_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned formal source exactly")
        object.__setattr__(self, "dependency_symbols", symbols)


@dataclass(frozen=True)
class ActualSignedInvariantCommonCurlAdmission:
    """Compose invariant-derived request jets with the actual common-curl theorem."""

    jets: ActualSignedInvariantJetsWitness
    common_curl: ActualSignedCommonCurlAdmission

    def __post_init__(self) -> None:
        if not isinstance(self.jets, ActualSignedInvariantJetsWitness):
            raise TypeError("jets must be an ActualSignedInvariantJetsWitness")
        if not isinstance(self.common_curl, ActualSignedCommonCurlAdmission):
            raise TypeError("common_curl must be an ActualSignedCommonCurlAdmission")
        failed = [name for name, ok in self.admission_checks().items() if not ok]
        if failed:
            raise ValueError(
                "uncertified invariant-to-common-curl composition: " + ", ".join(failed)
            )

    def admission_checks(self) -> dict[str, bool]:
        curl = self.common_curl.witness
        return {
            "pinned_repository": self.jets.lean_repository == PINNED_LEAN_REPOSITORY,
            "pinned_commit": self.jets.lean_commit == PINNED_LEAN_COMMIT,
            "pinned_jets_file": self.jets.lean_file == PINNED_JETS_FILE,
            "pinned_jets_theorem": self.jets.theorem == PINNED_JETS_THEOREM,
            "pinned_jets_dependencies": self.jets.dependency_symbols
            == REQUIRED_JETS_SYMBOLS,
            "theorem_evidence_not_sampled": self.jets.evidence_kind in _ACCEPTED_EVIDENCE,
            "cycle_analytic_invariant_certified": self.jets.cycle_analytic_invariant_certified
            is True,
            "full_request_identity_certified": self.jets.full_request_identity_certified
            is True,
            "geometry_strip_is_actual_strip_certified": self.jets.geometry_strip_is_actual_strip_certified
            is True,
            "phase_cell_instantiation_certified": self.jets.phase_cell_instantiation_certified
            is True,
            "beta_is_invariant_sigma_certified": self.jets.beta_is_invariant_sigma_certified
            is True,
            "jets_theorem_application_certified": self.jets.theorem_application_certified
            is True,
            "same_B": self.jets.B == curl.B,
            "same_N0": self.jets.N0 == curl.N0,
            "same_request": self.jets.request_repr == curl.request_repr,
            "same_beta": self.jets.beta_repr == curl.beta_repr,
            "common_curl_admitted": self.common_curl.formal_common_curl_theorem_admitted,
            "provenance_present": bool(self.jets.provenance),
        }

    @property
    def application_key(self) -> tuple[int, int, str, str, str, str]:
        """Stable identity for this composed formal application."""

        return (
            self.jets.B,
            self.jets.N0,
            self.jets.request_repr,
            self.jets.beta_repr,
            self.jets.application_id,
            self.common_curl.witness.application_id,
        )

    @property
    def invariant_supplies_common_curl_jets(self) -> bool:
        return all(self.admission_checks().values())

    @property
    def formal_divergence_zero_theorem_admitted(self) -> bool:
        return (
            self.invariant_supplies_common_curl_jets
            and self.common_curl.formal_divergence_zero_theorem_admitted
        )

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def full_request_jets_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_cycle_invariant_machine_verified(self) -> bool:
        return False

    @property
    def common_curl_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_signed_wave_values_materialized(self) -> bool:
        return False

    @property
    def paper_exact_divergence_free_wave_available(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
