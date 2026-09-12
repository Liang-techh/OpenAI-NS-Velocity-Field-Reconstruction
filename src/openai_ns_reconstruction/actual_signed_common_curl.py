"""Fail-closed theorem admission for the actual common signed curl wave.

The repository already contains coefficient-level curl algebra in
``curl_realization_algebra.py``.  The pinned formal source goes further:
``NavierStokes/ActualSignedCommonDynamics.lean`` proves
``common_curl_and_divergence`` for the actual signed-copy construction.  On a
point of ``fullStrip.domain`` and under the theorem's uniform-local-jet
hypothesis, the common corrected signed wave is realized as the cylindrical
curl of the common curl potential and its vector mode has zero cylindrical
divergence.

This module does not replay Lean and does not materialize any field.  It only
admits externally theorem-certified application metadata when the exact source,
declaration, hypotheses, actual-copy identities, and dependency symbols all
match the pinned formal construction.  Sampled, fitted, or numerical evidence
cannot pass this gate.

Successful admission therefore advances the formal Section 7 wave interface
without claiming that the repository has the paper's concrete amplitude/base
inputs, an actual theorem replay, or a paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral


PINNED_LEAN_REPOSITORY = "openai/NavierStokesAndEuler"
PINNED_LEAN_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_LEAN_FILE = "NavierStokes/ActualSignedCommonDynamics.lean"
PINNED_COMMON_CURL_THEOREM = (
    "ActualSignedCommonDynamics.common_curl_and_divergence"
)
PINNED_COPIES_DEFINITION = "ActualSignedOutputBounds.copies"
PINNED_PARAMETERS_DEFINITION = "ActualSignedStageControls.parameters"
PINNED_PHASE_CELL_DEFINITION = "ActualSignedStageControls.phaseCell"
PINNED_RAW_JETS_THEOREM = "ActualSignedCommonDynamics.rawJets"
PINNED_NATIVE_CURL_THEOREM = "ClosedNativeWaveIdentities.native_realizes_curl_at"
PINNED_NATIVE_DIVERGENCE_THEOREM = (
    "ClosedNativeWaveIdentities.native_divergence_zero_at"
)
PINNED_ZERO_GERM_THEOREM = (
    "LocalizedCurlRealization.native_identities_of_zero_germ"
)

REQUIRED_DEPENDENCY_SYMBOLS = (
    PINNED_COPIES_DEFINITION,
    PINNED_PARAMETERS_DEFINITION,
    PINNED_PHASE_CELL_DEFINITION,
    PINNED_RAW_JETS_THEOREM,
    PINNED_NATIVE_CURL_THEOREM,
    PINNED_NATIVE_DIVERGENCE_THEOREM,
    PINNED_ZERO_GERM_THEOREM,
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
class ActualSignedCommonCurlWitness:
    """Metadata for one application of the pinned common-curl theorem.

    ``signed_label_repr`` and ``point_repr`` are intentionally opaque.  Python
    must not invent an encoding for Lean's ``SignedLabel B N0`` or ``FullPoint``.
    Their identity is carried only so a later machine-linked export can bind the
    same application without changing this contract.
    """

    B: int
    N0: int
    harmonic: int
    signed_label_repr: str
    point_repr: str
    request_repr: str
    beta_repr: str
    application_id: str
    evidence_kind: str
    provenance: str
    uniform_local_jets_hypothesis_certified: bool
    full_strip_domain_membership_certified: bool
    actual_signed_copy_identity_certified: bool
    actual_parameter_identity_certified: bool
    theorem_application_certified: bool
    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_LEAN_FILE
    theorem: str = PINNED_COMMON_CURL_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_DEPENDENCY_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "B", _natural(self.B, "B"))
        object.__setattr__(self, "N0", _natural(self.N0, "N0"))
        object.__setattr__(self, "harmonic", _natural(self.harmonic, "harmonic"))
        for name in (
            "signed_label_repr",
            "point_repr",
            "request_repr",
            "beta_repr",
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
            "uniform_local_jets_hypothesis_certified",
            "full_strip_domain_membership_certified",
            "actual_signed_copy_identity_certified",
            "actual_parameter_identity_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)

        pins = (
            ("lean_repository", PINNED_LEAN_REPOSITORY),
            ("lean_commit", PINNED_LEAN_COMMIT),
            ("lean_file", PINNED_LEAN_FILE),
            ("theorem", PINNED_COMMON_CURL_THEOREM),
        )
        for name, expected in pins:
            object.__setattr__(self, name, _pinned(getattr(self, name), expected, name))

        if not isinstance(self.dependency_symbols, tuple):
            raise ValueError("dependency_symbols must be the pinned theorem-symbol tuple")
        symbols = tuple(
            _nonempty_text(symbol, "dependency symbol") for symbol in self.dependency_symbols
        )
        if symbols != REQUIRED_DEPENDENCY_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned formal source exactly")
        object.__setattr__(self, "dependency_symbols", symbols)

    @property
    def application_key(self) -> tuple[int, int, str, int, str, str]:
        """Stable opaque identity of this theorem application."""

        return (
            self.B,
            self.N0,
            self.signed_label_repr,
            self.harmonic,
            self.point_repr,
            self.application_id,
        )


@dataclass(frozen=True)
class ActualSignedCommonCurlAdmission:
    """Validated formal admission of one actual signed common-curl theorem use."""

    witness: ActualSignedCommonCurlWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, ActualSignedCommonCurlWitness):
            raise TypeError("witness must be an ActualSignedCommonCurlWitness")
        failed = [name for name, ok in self.admission_checks().items() if not ok]
        if failed:
            raise ValueError("uncertified actual signed common curl: " + ", ".join(failed))

    def admission_checks(self) -> dict[str, bool]:
        return {
            "pinned_repository": self.witness.lean_repository == PINNED_LEAN_REPOSITORY,
            "pinned_commit": self.witness.lean_commit == PINNED_LEAN_COMMIT,
            "pinned_file": self.witness.lean_file == PINNED_LEAN_FILE,
            "pinned_theorem": self.witness.theorem == PINNED_COMMON_CURL_THEOREM,
            "pinned_dependencies": self.witness.dependency_symbols
            == REQUIRED_DEPENDENCY_SYMBOLS,
            "theorem_evidence_not_sampled": self.witness.evidence_kind in _ACCEPTED_EVIDENCE,
            "uniform_local_jets_hypothesis_certified": self.witness.uniform_local_jets_hypothesis_certified
            is True,
            "full_strip_domain_membership_certified": self.witness.full_strip_domain_membership_certified
            is True,
            "actual_signed_copy_identity_certified": self.witness.actual_signed_copy_identity_certified
            is True,
            "actual_parameter_identity_certified": self.witness.actual_parameter_identity_certified
            is True,
            "theorem_application_certified": self.witness.theorem_application_certified is True,
            "provenance_present": bool(self.witness.provenance),
        }

    @property
    def application_key(self) -> tuple[int, int, str, int, str, str]:
        return self.witness.application_key

    @property
    def formal_common_curl_theorem_admitted(self) -> bool:
        return all(self.admission_checks().values())

    @property
    def formal_divergence_zero_theorem_admitted(self) -> bool:
        return all(self.admission_checks().values())

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def theorem_application_machine_replayed(self) -> bool:
        return False

    @property
    def actual_signed_wave_values_materialized(self) -> bool:
        return False

    @property
    def canonical_scope_link_verified(self) -> bool:
        return False

    @property
    def amplitude_ode_inputs_verified(self) -> bool:
        return False

    @property
    def paper_exact_divergence_free_wave_available(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
